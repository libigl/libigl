#include "frame_field.h"

#include <igl/tt.h>
#include <igl/edgetopology.h>
#include <igl/per_face_normals.h>
#include <igl/comiso/nrosy.h>
#include <iostream>

namespace igl
{

class FrameInterpolator
{
public:
  // Constraint type
  enum ConstraintType { FREE, SOFT, HARD };

  // Init
  IGL_INLINE FrameInterpolator(const Eigen::MatrixXd& _V, const Eigen::MatrixXi& _F);
  IGL_INLINE ~FrameInterpolator();

  // Generate the frame field
  IGL_INLINE void solve(bool hard_const);

  IGL_INLINE void setConstraint(const int fid, const Eigen::VectorXd& v, ConstraintType type = HARD);

  IGL_INLINE void setConstraintProportionalScale(const int fid, const double scale, ConstraintType type = HARD);

  // Set the ratio between smoothness and soft constraints (0 -> smoothness only, 1 -> soft constr only)
  IGL_INLINE void setSoftAlpha(double alpha);

  // Reset constraints (at least one constraint must be present or solve will fail)
  IGL_INLINE void resetConstraints();

  // Return the current field
  IGL_INLINE Eigen::MatrixXd getFieldPerFace();

  // Return the singularities
  IGL_INLINE Eigen::VectorXd getSingularityIndexPerVertex();


  // -------------- This is untested and unstable code

  IGL_INLINE void setConstraint_polar(const int fid, const Eigen::VectorXd& v, ConstraintType type);

  IGL_INLINE void interpolateSymmetric_polar(const bool bilaplacian);

  // Generate the frame field
  IGL_INLINE void solve_polar(const bool bilaplacian,bool hard_const);

  // Convert the frame field in the canonical representation
  IGL_INLINE void frame2canonical_polar(const Eigen::MatrixXd& TP, const Eigen::RowVectorXd& v, double& theta, Eigen::VectorXd& S);

  // Convert the canonical representation in a frame field
  IGL_INLINE void canonical2frame_polar(const Eigen::MatrixXd& TP, const double theta, const Eigen::VectorXd& S, Eigen::RowVectorXd& v);

  IGL_INLINE Eigen::MatrixXd getFieldPerFace_polar();

  IGL_INLINE void PolarDecomposition(Eigen::MatrixXd V, Eigen::MatrixXd& U, Eigen::MatrixXd& P);

  // Symmetric
  Eigen::MatrixXd S_polar;
  std::vector<ConstraintType> S_polar_c;

  // -------------------------------------------------

  // Face Topology
  Eigen::MatrixXi TT, TTi;

  // Two faces are consistent if their representative vector are taken modulo PI
  std::vector<bool> edge_consistency;
  Eigen::MatrixXi   edge_consistency_TT;

private:
  IGL_INLINE double mod2pi(double d);
  IGL_INLINE double modpi2(double d);
  IGL_INLINE double modpi(double d);

  // Convert a direction on the tangent space into an angle
  IGL_INLINE double vector2theta(const Eigen::MatrixXd& TP, const Eigen::RowVectorXd& v);

  // Convert an angle in a vector in the tangent space
  IGL_INLINE Eigen::RowVectorXd theta2vector(const Eigen::MatrixXd& TP, const double theta);

  // Convert the frame field in the canonical representation
  IGL_INLINE void frame2canonical(const Eigen::MatrixXd& TP, const Eigen::RowVectorXd& v, double& theta, double& beta, Eigen::RowVectorXd& d);

  // Convert the canonical representation in a frame field
  IGL_INLINE void canonical2frame(const Eigen::MatrixXd& TP, const double theta, const double beta, const Eigen::RowVectorXd& d, Eigen::RowVectorXd& v);

  // Interpolate the cross field (theta)
  IGL_INLINE void interpolateCross(bool hard_const);

  // Interpolate the skewness (beta)
  IGL_INLINE void interpolateSkewness();

  // Interpolate the scale (d)
  IGL_INLINE void interpolateScale();

  // Compute difference between reference frames
  IGL_INLINE void computek();

  // Compute edge consistency
  IGL_INLINE void compute_edge_consistency();

  // Cross field direction
  Eigen::VectorXd thetas;
  std::vector<ConstraintType> thetas_c;

  // Skewness
  Eigen::VectorXd betas;
  std::vector<ConstraintType> betas_c;

  // Scale
  Eigen::MatrixXd ds;
  std::vector<ConstraintType> ds_c;

  // Soft constraints strenght
  double softAlpha;

  // Edge Topology
  Eigen::MatrixXi EV, FE, EF;
  std::vector<bool> isBorderEdge;

  // Angle between two reference frames
  // R(k) * t0 = t1
  Eigen::VectorXd k;

  // Mesh
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  // Normals per face
  Eigen::MatrixXd N;

  // Singularity index
  Eigen::VectorXd singularityIndex;

  // Reference frame per triangle
  std::vector<Eigen::MatrixXd> TPs;

};




FrameInterpolator::FrameInterpolator(const Eigen::MatrixXd& _V, const Eigen::MatrixXi& _F)
{
  using namespace std;
  using namespace Eigen;

  V = _V;
  F = _F;

  assert(V.rows() > 0);
  assert(F.rows() > 0);


  // Generate topological relations
  igl::tt(V,F,TT,TTi);
  igl::edgetopology(V,F, EV, FE, EF);

  // Flag border edges
  isBorderEdge.resize(EV.rows());
  for(unsigned i=0; i<EV.rows(); ++i)
    isBorderEdge[i] = (EF(i,0) == -1) || ((EF(i,1) == -1));

  // Generate normals per face
  igl::per_face_normals(V, F, N);

  // Generate reference frames
  for(unsigned fid=0; fid<F.rows(); ++fid)
  {
    // First edge
    Vector3d e1 = V.row(F(fid,1)) - V.row(F(fid,0));
    e1.normalize();
    Vector3d e2 = N.row(fid);
    e2 = e2.cross(e1);
    e2.normalize();

    MatrixXd TP(2,3);
    TP << e1.transpose(), e2.transpose();
    TPs.push_back(TP);
  }

  // Reset the constraints
  resetConstraints();

  // Compute k, differences between reference frames
  computek();

  softAlpha = 0.5;

  // Alloc internal variables
  thetas            = VectorXd::Zero(F.rows());
//  thetas_c.resize(F.rows());

  betas            = VectorXd::Zero(F.rows());
//  betas_c.resize(F.rows());

  ds                = MatrixXd::Constant(F.rows(),2,1);
//  ds_c.resize(F.rows());

  S_polar = MatrixXd::Zero(F.rows(),3);
//  S_polar_c.resize(F.rows());

  singularityIndex = VectorXd::Zero(V.rows());

  compute_edge_consistency();

}

FrameInterpolator::~FrameInterpolator()
{

}

void FrameInterpolator::solve(bool hard_const)
{
    interpolateCross(hard_const);
    interpolateSkewness();
    interpolateScale();
}

void FrameInterpolator::setConstraint(const int fid, const Eigen::VectorXd& v, ConstraintType type)
{
  using namespace std;
  using namespace Eigen;
  
  double   t_;
  double   b_;
  RowVectorXd d_;
  frame2canonical(TPs[fid],v,t_,b_,d_);

//  cerr << "TP: " << endl << TPs[fid] << endl;
//  cerr << "v: " << endl << TPs[fid] << endl;
//  cerr << "t_: " << endl << TPs[fid] << endl;
//  cerr << "b_: " << endl << TPs[fid] << endl;
//  cerr << "d_: " << endl << TPs[fid] << endl;

  Eigen::RowVectorXd v2;
  canonical2frame(TPs[fid], t_, b_, d_, v2);
//  cerr << "Cosntraint: " << t_ << " " << b_ << " " << d_ << endl;
//  cerr << "Compare: " << endl << v.transpose() << endl << v2 << endl;

  thetas(fid)   = t_;
  thetas_c[fid] = type;

  betas(fid)    = b_;
  betas_c[fid]  = type;

  ds.row(fid)   = d_;
  ds_c[fid]     = type;
}

void FrameInterpolator::setConstraintProportionalScale(const int fid, const double scale, ConstraintType type)
{
  using namespace std;
  using namespace Eigen;
  

  assert(scale>0);
  ds.row(fid)   = ds.row(fid) * scale;
  ds_c[fid]     = type;
}

void FrameInterpolator::setSoftAlpha(double alpha)
{
  assert(alpha >= 0 && alpha <= 1);
  softAlpha = alpha;
}


Eigen::MatrixXd FrameInterpolator::getFieldPerFace()
{
  Eigen::MatrixXd R(F.rows(),6);
  for (unsigned i=0; i<F.rows(); ++i)
  {
    Eigen::RowVectorXd v;
    canonical2frame(TPs[i],thetas(i),betas(i),ds.row(i),v);
    R.row(i) = v;
  }
  return R;
}

Eigen::VectorXd FrameInterpolator::getSingularityIndexPerVertex()
{
  assert(0);
}


double FrameInterpolator::mod2pi(double d)
{
  while(d<0)
    d = d + (2.0*M_PI);

  return fmod(d, (2.0*M_PI));
}

double FrameInterpolator::modpi2(double d)
{
  while(d<0)
    d = d + (M_PI/2.0);

  return fmod(d, (M_PI/2.0));
}

double FrameInterpolator::modpi(double d)
{
  while(d<0)
    d = d + (M_PI);

  return fmod(d, (M_PI));
}


double FrameInterpolator::vector2theta(const Eigen::MatrixXd& TP, const Eigen::RowVectorXd& v)
{
  // Project onto the tangent plane
  Eigen::Vector2d vp = TP * v.transpose();

  // Convert to angle
  double theta = atan2(vp(1),vp(0));

//  cerr << v << endl << theta2vector(TP, theta) << endl;
//  assert((v.normalized() - theta2vector(TP, theta)).norm() < 10e-5);

  return theta;
}

Eigen::RowVectorXd FrameInterpolator::theta2vector(const Eigen::MatrixXd& TP, const double theta)
{
  Eigen::Vector2d vp(cos(theta),sin(theta));
  return vp.transpose() * TP;
}

void FrameInterpolator::frame2canonical(const Eigen::MatrixXd& TP, const Eigen::RowVectorXd& v, double& theta, double& beta, Eigen::RowVectorXd& d)
{
  using namespace std;
  using namespace Eigen;
  

  bool debug = false;
  if (debug)
  {
    cerr << "TP:" << endl << TP << endl;
    cerr << "v:" << endl << v << endl;
  }

  RowVectorXd v0 = v.segment<3>(0);
  RowVectorXd v1 = v.segment<3>(3);

  if (debug)
  {
    cerr << "v0:" << endl << v0 << endl;
    cerr << "v1:" << endl << v1 << endl;
  }

  // Find canonical cross field
  double theta0 = mod2pi(vector2theta(TP,v0));
  double theta1 = mod2pi(vector2theta(TP,v1));
  //cerr << mod2pi(theta0) << " " << mod2pi(theta1) << endl;

  // Find the closest pair
//  if (mod2pi(abs(theta0 + M_PI - theta1)) < mod2pi(abs(theta0 - theta1)))
//      theta0 = theta0 + M_PI;

  // Express theta1 as a displacement wrt theta0
  theta1 = mod2pi(theta1 - theta0);

  // Find the closest representative
  if (theta1 > M_PI)
    theta1 -= M_PI;
  assert(theta1 < M_PI);

  // Extract theta and beta
  theta = mod2pi(theta1/2.0 + theta0);
  beta  = theta1/2.0;

  //cerr << mod2pi(beta + theta) << " " << mod2pi(-beta + theta) << endl;

  d.resize(2);
  // Find the scaling factors
  d(0) = v0.norm();
  d(1) = v1.norm();

  if (debug)
  {
    cerr << "d:" << endl << d << endl;
    cerr << "thetavector:" << endl << theta2vector(TP, theta) << endl;
  }

}

void FrameInterpolator::canonical2frame(const Eigen::MatrixXd& TP, const double theta, const double beta, const Eigen::RowVectorXd& d, Eigen::RowVectorXd& v)
{
  using namespace std;
  using namespace Eigen;
  

  assert(d.cols() == 2 && d.rows() == 1);

  double theta0 = theta + M_PI - beta;
  double theta1 = theta + M_PI + beta;

  RowVectorXd v0 = theta2vector(TP,modpi(theta0)) * d(0);
  RowVectorXd v1 = theta2vector(TP,modpi(theta1)) * d(1);

//  cerr << "v0 : " << v0 << endl;
//  cerr << "v1 : " << v1 << endl;

  v.resize(6);
  v << v0, v1;
}

void FrameInterpolator::interpolateCross(bool hard_const)
{
  using namespace std;
  using namespace Eigen;

  NRosyField nrosy(V,F);
  nrosy.setSoftAlpha(softAlpha);

  for (unsigned i=0; i<F.rows(); ++i)
  {
    switch (thetas_c[i])
    {
      case FREE:
      break;
      case SOFT:
//      nrosy.setConstraintSoft(i,1,theta2vector(TPs[i],thetas(i)));
//      cerr << theta2vector(TPs[i],thetas(i)) << endl;
//      break;
      case HARD:
//      nrosy.setConstraintHard(i,theta2vector(TPs[i],thetas(i)));
        if (hard_const)
          nrosy.setConstraintHard(i,theta2vector(TPs[i],thetas(i)));
        else
          nrosy.setConstraintSoft(i,1,theta2vector(TPs[i],thetas(i)));
      break;
      default:
      assert(0);
    }
  }

  nrosy.solve(4);

  MatrixXd R = nrosy.getFieldPerFace();
//  std::cerr << "Cross: " << R << std::endl;

  assert(R.rows() == F.rows());
  for (unsigned i=0; i<F.rows(); ++i)
    if (thetas_c[i] != HARD)
      thetas(i) = vector2theta(TPs[i],R.row(i));
}

void FrameInterpolator::interpolateSkewness()
{
  using namespace std;
  using namespace Eigen;
  
  compute_edge_consistency();

  // Generate enriched Laplacian matrix (see notes)
  typedef Eigen::Triplet<double> triplet;
  std::vector<triplet> triplets;
  triplets.reserve(4*F.rows());

  VectorXd b = VectorXd::Zero(F.rows());

  // Build L and b
  for (unsigned eid=0; eid<EF.rows(); ++eid)
  {
    if (!isBorderEdge[eid])
    {
      int i = EF(eid,0);
      int j = EF(eid,1);
      double v = edge_consistency[eid]? 1 : -1;

      // Off-diagonal, symmetric
      triplets.push_back(triplet(i,j,v));
      triplets.push_back(triplet(j,i,v));
      // Diagonal
      triplets.push_back(triplet(i,i,-1));
      triplets.push_back(triplet(j,j,-1));

      if (!edge_consistency[eid])
      {
        b(i) -= M_PI/2.0; // CHECK
        b(j) -= M_PI/2.0; // CHECK
      }
    }
  }

  // Add soft constraints
  double w = 100000;
  for (unsigned fid=0; fid < F.rows(); ++fid)
  {
    if (betas_c[fid] != FREE)
    {
      triplets.push_back(triplet(fid,fid,w));
      b(fid) += w*betas(fid);
    }
  }

  SparseMatrix<double> L(F.rows(),F.rows());
  L.setFromTriplets(triplets.begin(), triplets.end());

  // Solve Lx = b;

  SimplicialLDLT<SparseMatrix<double> > solver;
  solver.compute(L);

  if(solver.info()!=Success)
  {
    std::cerr << "Cholesky failed - frame_interpolator.cpp" << std::endl;
    assert(0);
  }

  VectorXd x;
  x = solver.solve(b);

//  cerr << "Skewness: " << x << endl;

  if(solver.info()!=Success)
  {
    std::cerr << "Linear solve failed - frame_interpolator.cpp" << std::endl;
    assert(0);
  }

  // Copy back the result
  betas = x;
}

void FrameInterpolator::interpolateScale()
{
  using namespace std;
  using namespace Eigen;

  compute_edge_consistency();

  // Generate enriched Laplacian matrix (see notes)
  // the variables here are d1, d2
  typedef Eigen::Triplet<double> triplet;
  std::vector<triplet> triplets;
  triplets.reserve(4*F.rows()*2);

  VectorXd b = VectorXd::Zero(F.rows()*2);

  // Build L and b
  for (unsigned eid=0; eid<EF.rows(); ++eid)
  {
    if (!isBorderEdge[eid])
    {
      int d1  = EF(eid,0);
      int d2  = EF(eid,0)+F.rows();
      int d1p = EF(eid,1);
      int d2p = EF(eid,1)+F.rows();

      if (edge_consistency[eid])
      {
        triplets.push_back(triplet(d1 ,d1p,1));
        triplets.push_back(triplet(d2 ,d2p,1));
        triplets.push_back(triplet(d1p,d1 ,1));
        triplets.push_back(triplet(d2p,d2 ,1));
      }
      else
      {
        triplets.push_back(triplet(d1 ,d2p,1));
        triplets.push_back(triplet(d2 ,d1p,1));
        triplets.push_back(triplet(d1p,d2 ,1));
        triplets.push_back(triplet(d2p,d1 ,1));
      }

      // Diagonal
      triplets.push_back(triplet(d1,d1,  -1));
      triplets.push_back(triplet(d2,d2,  -1));
      triplets.push_back(triplet(d1p,d1p,-1));
      triplets.push_back(triplet(d2p,d2p,-1));
    }
  }

  // Add soft constraints
  double w = 100000;
  for (unsigned fid=0; fid < F.rows(); ++fid)
  {
    if (ds_c[fid] != FREE)
    {
      int d1 = fid;
      int d2 = fid + F.rows();

      triplets.push_back(triplet(d1,d1,w));
      b(d1) += w*ds(fid,0);

      triplets.push_back(triplet(d2,d2,w));
      b(d2) += w*ds(fid,1);
    }
  }

  SparseMatrix<double> L(F.rows()*2,F.rows()*2);
  L.setFromTriplets(triplets.begin(), triplets.end());

  // Solve Lx = b;

  SimplicialLDLT<SparseMatrix<double> > solver;
  solver.compute(L);

  if(solver.info()!=Success)
  {
    std::cerr << "Cholesky failed - frame_interpolator.cpp" << std::endl;
    assert(0);
  }

  VectorXd x;
  x = solver.solve(b);

//  cerr << "Scale: " << endl << x << endl;
  if(solver.info()!=Success)
  {
    std::cerr << "Linear solve failed - frame_interpolator.cpp" << std::endl;
    assert(0);
  }

  // Copy back the result
  ds << x.block(0, 0, ds.rows(), 1), x.block(ds.rows(), 0, ds.rows(), 1);
}

void FrameInterpolator::resetConstraints()
{
  thetas_c.resize(F.rows());
  betas_c.resize(F.rows());
  ds_c.resize(F.rows());
  S_polar_c.resize(F.rows());

  for(unsigned i=0; i<F.rows(); ++i)
  {
    thetas_c[i]  = FREE;
    betas_c[i]   = FREE;
    ds_c[i]      = FREE;
    S_polar_c[i] = FREE;
  }

}

void FrameInterpolator::compute_edge_consistency()
{
  using namespace std;
  using namespace Eigen;

  // Compute per-edge consistency
  edge_consistency.resize(EF.rows());
  edge_consistency_TT = MatrixXi::Constant(TT.rows(),3,-1);

  // For every non-border edge
  for (unsigned eid=0; eid<EF.rows(); ++eid)
  {
    if (!isBorderEdge[eid])
    {
      int fid0 = EF(eid,0);
      int fid1 = EF(eid,1);

      double theta0 = thetas(fid0);
      double theta1 = thetas(fid1);

      theta0 = theta0 + k(eid);

      double r = modpi(theta0-theta1);

      edge_consistency[eid] = r < M_PI/4.0 || r > 3*(M_PI/4.0);

      // Copy it into edge_consistency_TT
      int i1 = -1;
      int i2 = -1;
      for (unsigned i=0; i<3; ++i)
      {
        if (TT(fid0,i) == fid1)
          i1 = i;
        if (TT(fid1,i) == fid0)
          i2 = i;
      }
      assert(i1 != -1);
      assert(i2 != -1);

      edge_consistency_TT(fid0,i1) = edge_consistency[eid];
      edge_consistency_TT(fid1,i2) = edge_consistency[eid];
    }
  }
}

void FrameInterpolator::computek()
{
  using namespace std;
  using namespace Eigen;

  k.resize(EF.rows());

  // For every non-border edge
  for (unsigned eid=0; eid<EF.rows(); ++eid)
  {
    if (!isBorderEdge[eid])
    {
      int fid0 = EF(eid,0);
      int fid1 = EF(eid,1);

      Vector3d N0 = N.row(fid0);
      //Vector3d N1 = N.row(fid1);

      // find common edge on triangle 0 and 1
      int fid0_vc = -1;
      int fid1_vc = -1;
      for (unsigned i=0;i<3;++i)
      {
        if (EV(eid,0) == F(fid0,i))
          fid0_vc = i;
        if (EV(eid,1) == F(fid1,i))
          fid1_vc = i;
      }
      assert(fid0_vc != -1);
      assert(fid1_vc != -1);

      Vector3d common_edge = V.row(F(fid0,(fid0_vc+1)%3)) - V.row(F(fid0,fid0_vc));
      common_edge.normalize();

      // Map the two triangles in a new space where the common edge is the x axis and the N0 the z axis
      MatrixXd P(3,3);
      VectorXd o = V.row(F(fid0,fid0_vc));
      VectorXd tmp = -N0.cross(common_edge);
      P << common_edge, tmp, N0;
      P.transposeInPlace();


      MatrixXd V0(3,3);
      V0.row(0) = V.row(F(fid0,0)).transpose() -o;
      V0.row(1) = V.row(F(fid0,1)).transpose() -o;
      V0.row(2) = V.row(F(fid0,2)).transpose() -o;

      V0 = (P*V0.transpose()).transpose();

      assert(V0(0,2) < 10e-10);
      assert(V0(1,2) < 10e-10);
      assert(V0(2,2) < 10e-10);

      MatrixXd V1(3,3);
      V1.row(0) = V.row(F(fid1,0)).transpose() -o;
      V1.row(1) = V.row(F(fid1,1)).transpose() -o;
      V1.row(2) = V.row(F(fid1,2)).transpose() -o;
      V1 = (P*V1.transpose()).transpose();

      assert(V1(fid1_vc,2) < 10e-10);
      assert(V1((fid1_vc+1)%3,2) < 10e-10);

      // compute rotation R such that R * N1 = N0
      // i.e. map both triangles to the same plane
      double alpha = -atan2(V1((fid1_vc+2)%3,2),V1((fid1_vc+2)%3,1));

      MatrixXd R(3,3);
      R << 1,          0,            0,
           0, cos(alpha), -sin(alpha) ,
           0, sin(alpha),  cos(alpha);
      V1 = (R*V1.transpose()).transpose();

      assert(V1(0,2) < 10e-10);
      assert(V1(1,2) < 10e-10);
      assert(V1(2,2) < 10e-10);

      // measure the angle between the reference frames
      // k_ij is the angle between the triangle on the left and the one on the right
      VectorXd ref0 = V0.row(1) - V0.row(0);
      VectorXd ref1 = V1.row(1) - V1.row(0);

      ref0.normalize();
      ref1.normalize();

      double ktemp = atan2(ref1(1),ref1(0)) - atan2(ref0(1),ref0(0));

      // just to be sure, rotate ref0 using angle ktemp...
      MatrixXd R2(2,2);
      R2 << cos(ktemp), -sin(ktemp), sin(ktemp), cos(ktemp);

      tmp = R2*ref0.head<2>();

      assert(tmp(0) - ref1(0) < (0.000001));
      assert(tmp(1) - ref1(1) < (0.000001));

      k[eid] = ktemp;
    }
  }

}


  void FrameInterpolator::frame2canonical_polar(const Eigen::MatrixXd& TP, const Eigen::RowVectorXd& v, double& theta, Eigen::VectorXd& S_v)
{
  using namespace std;
  using namespace Eigen;

  RowVectorXd v0 = v.segment<3>(0);
  RowVectorXd v1 = v.segment<3>(3);

  // Project onto the tangent plane
  Vector2d vp0 = TP * v0.transpose();
  Vector2d vp1 = TP * v1.transpose();

  // Assemble matrix
  MatrixXd M(2,2);
  M << vp0, vp1;

  if (M.determinant() < 0)
    M.col(1) = -M.col(1);

  assert(M.determinant() > 0);

  // cerr << "M: " << M << endl;

  MatrixXd R,S;
  // PolarDecomposition
  PolarDecomposition(M,R,S);

//  cerr << M << endl;
//  cerr << "--------" << endl;
//  cerr << R*S << endl;

  // Finally, express the cross field as an angle
  theta = atan2(R(1,0),R(0,0)); // atan2(R(1,0),R(0,0));

  MatrixXd R2(2,2);
  R2 << cos(theta), -sin(theta), sin(theta), cos(theta);

//  cerr << "R: " << R << endl;
//  cerr << "Should be zero2: " << R2-R << endl;

  assert((R2-R).norm() < 10e-8);


  // Convert into rotation invariant form
  S = R * S * R.inverse();

  // cerr << "R" << endl << R << endl;
  // cerr << "S" << endl << S << endl;

  // Copy in vector form
  S_v = VectorXd(3);
  S_v << S(0,0), S(0,1), S(1,1);

}

  void FrameInterpolator::canonical2frame_polar(const Eigen::MatrixXd& TP, const double theta, const Eigen::VectorXd& S_v, Eigen::RowVectorXd& v)
{
  using namespace std;
  using namespace Eigen;

  assert(S_v.size() == 3);

  MatrixXd S(2,2);
  S << S_v(0), S_v(1), S_v(1), S_v(2);

  // Convert angle in vector in the tangent plane
  // Vector2d vp(cos(theta),sin(theta));

  // First reconstruct R
  MatrixXd R(2,2);
  // R.col(0) << cos(theta), sin(theta);
  // R(0,1) = -R(1,0);
  // R(1,1) =  R(0,0);
  R << cos(theta), -sin(theta), sin(theta), cos(theta);


  // cerr << "R2" << endl << R << endl;
  // cerr << "S2" << endl << S << endl;

  // Multiply to reconstruct
  //MatrixXd M = R * S;

  // Rotation invariant reconstruction
  MatrixXd M = S * R;

  // cerr << "M2: " << M << endl;

  Vector2d vp0(M(0,0),M(1,0));
  Vector2d vp1(M(0,1),M(1,1));

  // Unproject the vectors
  RowVectorXd v0 = vp0.transpose() * TP;
  RowVectorXd v1 = vp1.transpose() * TP;

  v.resize(6);
  v << v0, v1;
}

void FrameInterpolator::solve_polar(const bool bilaplacian, bool hard_const)
{
    interpolateCross(hard_const);
    //interpolateSymmetric_polar(bilaplacian);
  interpolateSymmetric_polar(false);
}

void FrameInterpolator::interpolateSymmetric_polar(const bool bilaplacian)
{
  using namespace std;
  using namespace Eigen;

  //compute_edge_consistency();

  // Generate uniform Laplacian matrix (see notes)
  typedef Eigen::Triplet<double> triplet;
  std::vector<triplet> triplets;

  // Variables are stacked as x1,y1,z1,x2,y2,z2
  triplets.reserve(3*4*F.rows());

  MatrixXd b = MatrixXd::Zero(3*F.rows(),1);

  // Build L and b
  for (unsigned eid=0; eid<EF.rows(); ++eid)
  {
    if (!isBorderEdge[eid])
    {
      for (int z=0;z<2;++z)
      {
        // W = [w_a, w_b
        //      w_b, w_c]
        //

        // It is not symmetric
        int i    = EF(eid,z==0?0:1);
        int j    = EF(eid,z==0?1:0);

        int w_a_0 = (i*3)+0;
        int w_b_0 = (i*3)+1;
        int w_c_0 = (i*3)+2;

        int w_a_1 = (j*3)+0;
        int w_b_1 = (j*3)+1;
        int w_c_1 = (j*3)+2;

        // Rotation to change frame
        double r_a =  cos(z==1?k(eid):-k(eid));
        double r_b = -sin(z==1?k(eid):-k(eid));
        double r_c =  sin(z==1?k(eid):-k(eid));
        double r_d =  cos(z==1?k(eid):-k(eid));

        if (true)
        {
          if (true)
          {
            // First term
            // w_a_0 = r_a^2 w_a_1 + 2 r_a r_b w_b_1 + r_b^2 w_c_1 = 0
            triplets.push_back(triplet(w_a_0,w_a_0,                -1 ));
            triplets.push_back(triplet(w_a_0,w_a_1,           r_a*r_a ));
            triplets.push_back(triplet(w_a_0,w_b_1,       2 * r_a*r_b ));
            triplets.push_back(triplet(w_a_0,w_c_1,           r_b*r_b ));

            // Second term
            // w_b_0 = r_a r_c w_a + (r_b r_c + r_a r_d) w_b + r_b r_d w_c
            triplets.push_back(triplet(w_b_0,w_b_0,                -1 ));
            triplets.push_back(triplet(w_b_0,w_a_1,           r_a*r_c ));
            triplets.push_back(triplet(w_b_0,w_b_1, r_b*r_c + r_a*r_d ));
            triplets.push_back(triplet(w_b_0,w_c_1,           r_b*r_d ));

            // Third term
            // w_c_0 = r_c^2 w_a + 2 r_c r_d w_b +  r_d^2 w_c
            triplets.push_back(triplet(w_c_0,w_c_0,                -1 ));
            triplets.push_back(triplet(w_c_0,w_a_1,           r_c*r_c ));
            triplets.push_back(triplet(w_c_0,w_b_1,       2 * r_c*r_d ));
            triplets.push_back(triplet(w_c_0,w_c_1,           r_d*r_d ));
          }
          else
          {
            // First term
            // w_a_0 = r_a^2 w_a_1 + 2 r_a r_b w_b_1 + r_b^2 w_c_1 = 0
            triplets.push_back(triplet(w_a_0,w_a_1,           r_a*r_a ));
            triplets.push_back(triplet(w_a_0,w_b_1,       2 * r_a*r_b ));
            triplets.push_back(triplet(w_a_0,w_c_1,           r_b*r_b ));
            triplets.push_back(triplet(w_a_0,w_a_0, -          r_a*r_a ));
            triplets.push_back(triplet(w_a_0,w_a_0, -      2 * r_a*r_b ));
            triplets.push_back(triplet(w_a_0,w_a_0, -          r_b*r_b ));

            // Second term
            // w_b_0 = r_a r_c w_a + (r_b r_c + r_a r_d) w_b + r_b r_d w_c
            triplets.push_back(triplet(w_b_0,w_a_1,           r_a*r_c ));
            triplets.push_back(triplet(w_b_0,w_b_1, r_b*r_c + r_a*r_d ));
            triplets.push_back(triplet(w_b_0,w_c_1,           r_b*r_d ));
            triplets.push_back(triplet(w_b_0,w_b_0, -          r_a*r_c ));
            triplets.push_back(triplet(w_b_0,w_b_0, -r_b*r_c + r_a*r_d ));
            triplets.push_back(triplet(w_b_0,w_b_0, -          r_b*r_d ));

            // Third term
            // w_c_0 = r_c^2 w_a + 2 r_c r_d w_b +  r_d^2 w_c
            triplets.push_back(triplet(w_c_0,w_a_1,           r_c*r_c ));
            triplets.push_back(triplet(w_c_0,w_b_1,       2 * r_c*r_d ));
            triplets.push_back(triplet(w_c_0,w_c_1,           r_d*r_d ));
            triplets.push_back(triplet(w_c_0,w_c_0, -          r_c*r_c ));
            triplets.push_back(triplet(w_c_0,w_c_0, -      2 * r_c*r_d ));
            triplets.push_back(triplet(w_c_0,w_c_0, -          r_d*r_d ));

          }
        }
        else
        {
          // Simple Laplacian
          // w_a_0 = r_a^2 w_a_1 + 2 r_a r_b w_b_1 + r_b^2 w_c_1 = 0
          triplets.push_back(triplet(w_a_0,w_a_0,                -1 ));
          triplets.push_back(triplet(w_a_0,w_a_1,                 1 ));

          // Second term
          // w_b_0 = r_a r_c w_a + (r_b r_c + r_a r_d) w_b + r_b r_d w_c
          triplets.push_back(triplet(w_b_0,w_b_0,                -1 ));
          triplets.push_back(triplet(w_b_0,w_b_1,                 1 ));

          // Third term
          // w_c_0 = r_c^2 w_a + 2 r_c r_d w_b +  r_d^2 w_c
          triplets.push_back(triplet(w_c_0,w_c_0,                -1 ));
          triplets.push_back(triplet(w_c_0,w_c_1,                 1 ));
        }
      }
    }
  }

  SparseMatrix<double> L(3*F.rows(),3*F.rows());
  L.setFromTriplets(triplets.begin(), triplets.end());

  triplets.clear();

  // Add soft constraints
  double w = 100000;
  for (unsigned fid=0; fid < F.rows(); ++fid)
  {
    if (S_polar_c[fid] != FREE)
    {
      for (unsigned i=0;i<3;++i)
      {
        triplets.push_back(triplet(3*fid + i,3*fid + i,w));
        b(3*fid + i) += w*S_polar(fid,i);
      }
    }
  }

  SparseMatrix<double> soft(3*F.rows(),3*F.rows());
  soft.setFromTriplets(triplets.begin(), triplets.end());

  SparseMatrix<double> M;

  if (!bilaplacian)
    M = L + soft;
  else
    M = L * L + soft;

  // Solve Lx = b;

//  MatrixXd Dense(L);
//
//  cerr << "Is sym:" << (Dense - Dense.transpose()).maxCoeff() << endl;

//  SimplicialLDLT<SparseMatrix<double> > solver;
  SparseLU<SparseMatrix<double> > solver;

  solver.compute(M);

  if(solver.info()!=Success)
  {
    std::cerr << "Cholesky failed - frame_interpolator.cpp" << std::endl;
    assert(0);
  }

  MatrixXd x;
  x = solver.solve(b);

//  cerr << "Skewness: " << x << endl;

  if(solver.info()!=Success)
  {
    std::cerr << "Linear solve failed - frame_interpolator.cpp" << std::endl;
    assert(0);
  }

  S_polar = MatrixXd::Zero(F.rows(),3);

  // Copy back the result
  for (unsigned i=0;i<F.rows();++i)
    S_polar.row(i) << x(i*3+0), x(i*3+1), x(i*3+2);

  for (unsigned i=0;i<F.rows();++i)
  {
    if (false)
    {
      MatrixXd T_POLAR(2,2);
      T_POLAR << x(i*3+0), x(i*3+1), x(i*3+1), x(i*3+2);

      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> S_es(T_POLAR);
      VectorXd eivals = S_es.eigenvalues();

//      cerr << "eivals:" << eivals << endl;
      double t = eivals.minCoeff();
      if (t < 0)
      {
        cerr << "NOOOOOOOO!" << endl;
        exit(1);
      }
//      cerr << S_polar << endl;
    }
  }
}

void FrameInterpolator::setConstraint_polar(const int fid, const Eigen::VectorXd& v, ConstraintType type)
{
  using namespace std;
  using namespace Eigen;

  double   t_;
  VectorXd S_;

  frame2canonical_polar(TPs[fid],v,t_,S_);

  Eigen::RowVectorXd v2;
  canonical2frame_polar(TPs[fid], t_, S_, v2);


//  cerr << "Constraint: " << t_ << " " << S_ << endl;
//  cerr << "Compare:    " << endl << v.transpose() << endl << v2 << endl;
//  cerr << "Diff:    " << endl << v.transpose() - v2 << endl;
//  cerr << "Diff:    " << endl << v.transpose() + v2 << endl;

  thetas(fid)   = t_;
  thetas_c[fid] = type;

  S_polar.row(fid) = S_;
  S_polar_c[fid]   = type;

}

Eigen::MatrixXd FrameInterpolator::getFieldPerFace_polar()
{
  using namespace std;
  using namespace Eigen;

  MatrixXd R(F.rows(),6);
  for (unsigned i=0; i<F.rows(); ++i)
  {
    RowVectorXd v;
    canonical2frame_polar(TPs[i],thetas(i),S_polar.row(i),v);
    R.row(i) = v;
  }
  return R;
}

  void FrameInterpolator::PolarDecomposition(Eigen::MatrixXd V, Eigen::MatrixXd& U, Eigen::MatrixXd& P)
{
  using namespace std;
  using namespace Eigen;

  // Polar Decomposition
  JacobiSVD<MatrixXd> svd(V,Eigen::ComputeFullU | Eigen::ComputeFullV);

  U = svd.matrixU() * svd.matrixV().transpose();
  P = svd.matrixV() * svd.singularValues().asDiagonal() * svd.matrixV().transpose();

  // If not SPD, change sign of both
  // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> S_es(P);
  // if (S_es.eigenvalues().minCoeff() < 0)
  // {
  //   U = -U;
  //   P = -P;
  // }
}

}


IGL_INLINE void igl::frame_field(
                                 const Eigen::MatrixXd& V,
                                 const Eigen::MatrixXi& F,
                                 const Eigen::VectorXi& b,
                                 const Eigen::MatrixXd& bc1,
                                 const Eigen::MatrixXd& bc2,
                                 Eigen::MatrixXd& F1,
                                 Eigen::MatrixXd& F2,
                                 Eigen::VectorXd& S
                                 )

{
  using namespace std;
  using namespace Eigen;

  assert(b.size() > 0);
  
  // Init Solver
  FrameInterpolator field(V,F);
  
  for (unsigned i=0; i<b.size(); ++i)
  {
    VectorXd t(6); t << bc1.row(i).transpose(), bc2.row(i).transpose();
    field.setConstraint_polar(b(i), t,FrameInterpolator::SOFT);
  }
  
  // Solve
  field.solve_polar(false,true);
  
  // Copy back
  MatrixXd R = field.getFieldPerFace_polar();
  F1 = R.block(0, 0, R.rows(), 3);
  F2 = R.block(0, 3, R.rows(), 3);
}

