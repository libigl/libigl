// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include <igl/angle_bound_frame_fields.h>
#include <igl/edge_topology.h>
#include <igl/local_basis.h>
#include <igl/sparse.h>
#include <igl/speye.h>
#include <igl/slice.h>
#include <igl/polyroots.h>
#include <igl/colon.h>
#include <Eigen/Sparse>

#include <iostream>

namespace igl {

  template <typename DerivedV, typename DerivedF>
  class AngleBoundFFSolverData
  {
    public:
      const Eigen::PlainObjectBase<DerivedV> &V; int numV;
      const Eigen::PlainObjectBase<DerivedF> &F; int numF;

      Eigen::MatrixXi EV; int numE;
      Eigen::MatrixXi F2E;
      Eigen::MatrixXi E2F;
      Eigen::VectorXd K;

      Eigen::VectorXi isBorderEdge;
      int numInteriorEdges;
      Eigen::Matrix<int,Eigen::Dynamic,2> E2F_int;
      Eigen::VectorXi indInteriorToFull;
      Eigen::VectorXi indFullToInterior;

      DerivedV B1, B2, FN;

      //laplacians
      Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar>> DDA, DDB;

  private:
    IGL_INLINE void computeLaplacians();
    IGL_INLINE void computek();
    IGL_INLINE void computeCoefficientLaplacian(int n, Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > &D);
    IGL_INLINE void precomputeInteriorEdges();

public:
      IGL_INLINE AngleBoundFFSolverData(const Eigen::PlainObjectBase<DerivedV> &_V,
                                   const Eigen::PlainObjectBase<DerivedF> &_F);
  };

  template <typename DerivedV, typename DerivedF, typename DerivedO>
  class AngleBoundFFSolver
  {
  public:
    IGL_INLINE AngleBoundFFSolver(const AngleBoundFFSolverData<DerivedV, DerivedF> &_data,
                                  const typename DerivedV::Scalar &_thetaMin = 30,
                                 int _maxIter = 50,
                                 const typename DerivedV::Scalar &_lambdaInit = 100,
                                 const typename DerivedV::Scalar &_lambdaMultFactor = 1.01,
                                const bool _doHardConstraints = false);
    IGL_INLINE bool solve(const Eigen::VectorXi &isConstrained,
                          const Eigen::PlainObjectBase<DerivedO> &initialSolution,
                          Eigen::PlainObjectBase<DerivedO> &output,
                          typename DerivedV::Scalar *lambdaOut = NULL);

  private:

    const AngleBoundFFSolverData<DerivedV, DerivedF> &data;

    //polyVF data
    Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> Acoeff, Bcoeff;
    Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 2> pvU, pvV;
    typename DerivedV::Scalar lambda;

    //parameters
    typename DerivedV::Scalar lambdaInit,lambdaMultFactor;
    int maxIter;
    typename DerivedV::Scalar thetaMin;
    bool doHardConstraints;

    typename DerivedV::Scalar computeAngle(const std::complex<typename DerivedV::Scalar> &u,
                                           const std::complex<typename DerivedV::Scalar> &v);
//    IGL_INLINE void computeAngles(Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 1> &angles);

    IGL_INLINE int getNumOutOfBounds();

    IGL_INLINE void rotateAroundBisector(const std::complex<typename DerivedV::Scalar> &uin,
                         const std::complex<typename DerivedV::Scalar> &vin,
                         const typename DerivedV::Scalar theta,
                         std::complex<typename DerivedV::Scalar> &uout,
                         std::complex<typename DerivedV::Scalar> &vout);

    IGL_INLINE void localStep();

    IGL_INLINE void globalStep(const Eigen::Matrix<int, Eigen::Dynamic, 1>  &isConstrained,
                               const Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1>  &Ak,
                               const Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1>  &Bk);

    IGL_INLINE void minQuadWithKnownMini(const Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > &Q,
                         const Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > &f,
                         const Eigen::VectorXi isConstrained,
                         const Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> &xknown,
                                         Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> &x);
    IGL_INLINE void setFieldFromCoefficients();
    IGL_INLINE void setCoefficientsFromField();

  };
}

//Implementation
/***************************** Data ***********************************/

template <typename DerivedV, typename DerivedF>
IGL_INLINE igl::AngleBoundFFSolverData<DerivedV, DerivedF>::
AngleBoundFFSolverData(const Eigen::PlainObjectBase<DerivedV> &_V,
                  const Eigen::PlainObjectBase<DerivedF> &_F):
V(_V),
numV(_V.rows()),
F(_F),
numF(_F.rows())
{
  igl::edge_topology(V,F,EV,F2E,E2F);
  numE = EV.rows();

  precomputeInteriorEdges();

  igl::local_basis(V,F,B1,B2,FN);

  computek();

  computeLaplacians();

};


template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::AngleBoundFFSolverData<DerivedV, DerivedF>::computeLaplacians()
{
  computeCoefficientLaplacian(2, DDA);

  computeCoefficientLaplacian(4, DDB);
}

template<typename DerivedV, typename DerivedF>
IGL_INLINE void igl::AngleBoundFFSolverData<DerivedV, DerivedF>::
precomputeInteriorEdges()
{
  // Flag border edges
  numInteriorEdges = 0;
  isBorderEdge.setZero(numE,1);
  indFullToInterior = Eigen::VectorXi::Constant(numE,-1);

  for(unsigned i=0; i<numE; ++i)
  {
    if ((E2F(i,0) == -1) || ((E2F(i,1) == -1)))
      isBorderEdge[i] = 1;
    else
    {
      indFullToInterior[i] = numInteriorEdges;
      numInteriorEdges++;
    }
  }

  E2F_int.resize(numInteriorEdges, 2);
  indInteriorToFull.setZero(numInteriorEdges,1);
  int ii = 0;
  for (int k=0; k<numE; ++k)
  {
    if (isBorderEdge[k])
      continue;
    E2F_int.row(ii) = E2F.row(k);
    indInteriorToFull[ii] = k;
    ii++;
  }

}



template<typename DerivedV, typename DerivedF>
IGL_INLINE void igl::AngleBoundFFSolverData<DerivedV, DerivedF>::
computeCoefficientLaplacian(int n, Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > &D)
{
  std::vector<Eigen::Triplet<std::complex<typename DerivedV::Scalar> >> tripletList;

  // For every non-border edge
  for (unsigned eid=0; eid<numE; ++eid)
  {
    if (!isBorderEdge[eid])
    {
      int fid0 = E2F(eid,0);
      int fid1 = E2F(eid,1);

      tripletList.push_back(Eigen::Triplet<std::complex<typename DerivedV::Scalar> >(fid0,
                                                                                     fid0,
                                                                                     std::complex<typename DerivedV::Scalar>(1.)));
      tripletList.push_back(Eigen::Triplet<std::complex<typename DerivedV::Scalar> >(fid1,
                                                                                     fid1,
                                                                                     std::complex<typename DerivedV::Scalar>(1.)));
      tripletList.push_back(Eigen::Triplet<std::complex<typename DerivedV::Scalar> >(fid0,
                                                                                     fid1,
                                                                                     -1.*std::polar(1.,-1.*n*K[eid])));
      tripletList.push_back(Eigen::Triplet<std::complex<typename DerivedV::Scalar> >(fid1,
                                                                                     fid0,
                                                                                     -1.*std::polar(1.,1.*n*K[eid])));

    }
  }
  D.resize(numF,numF);
  D.setFromTriplets(tripletList.begin(), tripletList.end());


}

template<typename DerivedV, typename DerivedF>
IGL_INLINE void igl::AngleBoundFFSolverData<DerivedV, DerivedF>::
computek()
{
  K.setZero(numE);
  // For every non-border edge
  for (unsigned eid=0; eid<numE; ++eid)
  {
    if (!isBorderEdge[eid])
    {
      int fid0 = E2F(eid,0);
      int fid1 = E2F(eid,1);

      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> N0 = FN.row(fid0);
      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> N1 = FN.row(fid1);

      // find common edge on triangle 0 and 1
      int fid0_vc = -1;
      int fid1_vc = -1;
      for (unsigned i=0;i<3;++i)
      {
        if (F2E(fid0,i) == eid)
          fid0_vc = i;
        if (F2E(fid1,i) == eid)
          fid1_vc = i;
      }
      assert(fid0_vc != -1);
      assert(fid1_vc != -1);

      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> common_edge = V.row(F(fid0,(fid0_vc+1)%3)) - V.row(F(fid0,fid0_vc));
      common_edge.normalize();

      // Map the two triangles in a new space where the common edge is the x axis and the N0 the z axis
      Eigen::Matrix<typename DerivedV::Scalar, 3, 3> P;
      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> o = V.row(F(fid0,fid0_vc));
      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> tmp = -N0.cross(common_edge);
      P << common_edge, tmp, N0;
      //      P.transposeInPlace();


      Eigen::Matrix<typename DerivedV::Scalar, 3, 3> V0;
      V0.row(0) = V.row(F(fid0,0)) -o;
      V0.row(1) = V.row(F(fid0,1)) -o;
      V0.row(2) = V.row(F(fid0,2)) -o;

      V0 = (P*V0.transpose()).transpose();

      Eigen::Matrix<typename DerivedV::Scalar, 3, 3> V1;
      V1.row(0) = V.row(F(fid1,0)) -o;
      V1.row(1) = V.row(F(fid1,1)) -o;
      V1.row(2) = V.row(F(fid1,2)) -o;
      V1 = (P*V1.transpose()).transpose();

      // compute rotation R such that R * N1 = N0
      // i.e. map both triangles to the same plane
      double alpha = -atan2(V1((fid1_vc+2)%3,2),V1((fid1_vc+2)%3,1));

      Eigen::Matrix<typename DerivedV::Scalar, 3, 3> R;
      R << 1,          0,            0,
      0, cos(alpha), -sin(alpha) ,
      0, sin(alpha),  cos(alpha);
      V1 = (R*V1.transpose()).transpose();

      // measure the angle between the reference frames
      // k_ij is the angle between the triangle on the left and the one on the right
      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> ref0 = V0.row(1) - V0.row(0);
      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> ref1 = V1.row(1) - V1.row(0);

      ref0.normalize();
      ref1.normalize();

      double ktemp = atan2(ref1(1),ref1(0)) - atan2(ref0(1),ref0(0));

      // just to be sure, rotate ref0 using angle ktemp...
      Eigen::Matrix<typename DerivedV::Scalar, 2, 2> R2;
      R2 << cos(ktemp), -sin(ktemp), sin(ktemp), cos(ktemp);

      Eigen::Matrix<typename DerivedV::Scalar, 1, 2> tmp1 = R2*(ref0.head(2)).transpose();

      K[eid] = ktemp;
    }
  }

}


/***************************** Solver ***********************************/
template <typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE igl::AngleBoundFFSolver<DerivedV, DerivedF, DerivedO>::
AngleBoundFFSolver(const AngleBoundFFSolverData<DerivedV, DerivedF> &_data,
                   const typename DerivedV::Scalar &_thetaMin,
                  int _maxIter,
                  const typename DerivedV::Scalar &_lambdaInit,
                  const typename DerivedV::Scalar &_lambdaMultFactor,
                   const bool _doHardConstraints):
data(_data),
lambdaInit(_lambdaInit),
maxIter(_maxIter),
lambdaMultFactor(_lambdaMultFactor),
doHardConstraints(_doHardConstraints),
thetaMin(_thetaMin)
{
  Acoeff.resize(data.numF,1);
  Bcoeff.resize(data.numF,1);
  pvU.setZero(data.numF, 2);
  pvV.setZero(data.numF, 2);
};

template<typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE void igl::AngleBoundFFSolver<DerivedV, DerivedF, DerivedO>::
rotateAroundBisector(const std::complex<typename DerivedV::Scalar> &uin,
                          const std::complex<typename DerivedV::Scalar> &vin,
                          const typename DerivedV::Scalar diff,
                          std::complex<typename DerivedV::Scalar> &uout,
                          std::complex<typename DerivedV::Scalar> &vout)
{
  //rotate 2D complex vectors u and v around their bisector so that their
  //angle is at least theta

  uout = uin;
  vout = vin;
  typename DerivedV::Scalar au = arg(uin);
  typename DerivedV::Scalar av = arg(vin);
  if (au<av)
  {
    uout = std::polar (1.0,-.5*diff)*uin;
    vout = std::polar (1.0, .5*diff)*vin;
  }
  else
  {
    uout = std::polar (1.0, .5*diff)*uin;
    vout = std::polar (1.0,-.5*diff)*vin;
  }

}


template<typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE void igl::AngleBoundFFSolver<DerivedV, DerivedF, DerivedO>::
localStep()
{
  for (int j =0; j<data.numF; ++j)
  {

    std::complex<typename DerivedV::Scalar> u(pvU(j,0),pvU(j,1));
    std::complex<typename DerivedV::Scalar> v(pvV(j,0),pvV(j,1));

    typename DerivedV::Scalar current_angle = computeAngle(u, v);
    if (current_angle<thetaMin*M_PI/180)
    {
      // bring all to 1st or 4th quarter plane
      if ((arg(u)>=0.5*M_PI || arg(u)<-0.5*M_PI ))
        u = -u;
      if ((arg(v)>=0.5*M_PI || arg(v)<-0.5*M_PI ))
        v = -v;
      assert(fabs(computeAngle(u, v) - current_angle)<1e-5);

      if ( fabs(arg(u) - arg(v)) >0.5*M_PI )
        v = -v;
      assert(fabs(computeAngle(u, v) - current_angle)<1e-5);

      std::complex<typename DerivedV::Scalar> u1, v1;
      typename DerivedV::Scalar diff = thetaMin*M_PI/180 - current_angle + 1e-6;
      rotateAroundBisector(u, v, diff, u1, v1);

//      if (computeAngle(u1, v1)<thetaMin*M_PI/180)
//      {
//        std::cerr<<"u = ["<<real(u)<<","<<imag(u)<< "]; v= ["<<real(v)<<","<<imag(v)<<"];"<<std::endl;
//        std::cerr<<"u1 = ["<<real(u1)<<","<<imag(u1)<< "]; v1= ["<<real(v1)<<","<<imag(v1)<<"];"<<std::endl;
//        std::cerr<<"current_angle = "<<current_angle<<std::endl;
//        std::cerr<<"aout = "<<computeAngle(u1, v1)<< "; theta= "<<thetaMin*M_PI/180<<";"<<std::endl;
//      }
//      assert(computeAngle(u1, v1)>=thetaMin*M_PI/180);


      pvU.row(j) << real(u1),imag(u1);
      pvV.row(j) << real(v1),imag(v1);
    }
  }

}


//
//template<typename DerivedV, typename DerivedF, typename DerivedO>
//IGL_INLINE void igl::AngleBoundFFSolver<DerivedV, DerivedF, DerivedO>::
//computeAngles(Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 1> &angles)
//{
//  angles.resize(data.numF,1);
//  for (int i =0; i<data.numF; ++i)
//  {
//    std::complex<typename DerivedV::Scalar> u(pvU(i,0),pvU(i,1));
//    std::complex<typename DerivedV::Scalar> v(pvV(i,0),pvV(i,1));
//    angles[i] = fabs(arg(u) - arg(v));
//    if (angles[i]>M_PI)
//      angles[i] = 2*M_PI-angles[i];
//    if (angles[i]>.5*M_PI)
//      angles[i] = M_PI-angles[i];
//  }
//}

template<typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE typename DerivedV::Scalar igl::AngleBoundFFSolver<DerivedV, DerivedF, DerivedO>::
computeAngle(const std::complex<typename DerivedV::Scalar> &u,
             const std::complex<typename DerivedV::Scalar> &v)
{
  typename DerivedV::Scalar angle = std::min(fabs(arg(u*conj(v))), fabs(arg(u*conj(-v))));

//  typename DerivedV::Scalar angle;
//  typename DerivedV::Scalar a1 = fabs(arg(u*conj(v)));
//  typename DerivedV::Scalar a2 = fabs(arg(u*conj(-v)));
//  if (a1 < a2)
//    angle = a1;
//  else
//  {
//    angle = a2; v = -v;
//  }

//  typename DerivedV::Scalar angle = fabs(arg(u) - arg(v));
//  if (angle>M_PI)
//  {
//    u = -u;
//    angle = fabs(arg(u) - arg(v));
//  };
//
//  if (angle>.5*M_PI)
//  {
//    v = -v;
//    angle = fabs(arg(u) - arg(v));
//  };
//
//  assert(fabs(angle-angle1)<1e-6);

//  if (angle>M_PI)
//    angle = 2*M_PI-angle;
//  if (angle>.5*M_PI)
//    angle = M_PI-angle;

//  typename DerivedV::Scalar angle = fabs(arg(u) - arg(v));
//    if (angle>M_PI)
//      angle = 2*M_PI-angle;
//    if (angle>.5*M_PI)
//      angle = M_PI-angle;

  assert(angle <= .5*M_PI && angle >0);

  return angle;
}


template<typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE int igl::AngleBoundFFSolver<DerivedV, DerivedF, DerivedO>::
getNumOutOfBounds()
{
  Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 1> angles;
//  computeAngles(angles);
  int numOoB = 0;
  for (int i =0; i<data.numF; ++i)
  {
    std::complex<typename DerivedV::Scalar> u(pvU(i,0),pvU(i,1));
    std::complex<typename DerivedV::Scalar> v(pvV(i,0),pvV(i,1));
    typename DerivedV::Scalar angle = computeAngle(u,v);
//    if (angles[i] <thetaMin*M_PI/180)
    if (angle <thetaMin*M_PI/180)
      numOoB ++;
  }
  return numOoB;
}

template<typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE void igl::AngleBoundFFSolver<DerivedV, DerivedF, DerivedO>::
setCoefficientsFromField()
{
  for (int i = 0; i <data.numF; ++i)
  {
    std::complex<typename DerivedV::Scalar> u(pvU(i,0),pvU(i,1));
    std::complex<typename DerivedV::Scalar> v(pvV(i,0),pvV(i,1));
    Acoeff(i) = u*u+v*v;
    Bcoeff(i) = u*u*v*v;
  }
}


template<typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE void igl::AngleBoundFFSolver<DerivedV, DerivedF, DerivedO>::
globalStep(const Eigen::Matrix<int, Eigen::Dynamic, 1>  &isConstrained,
           const Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1>  &Ak,
           const Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1>  &Bk)
{
  setCoefficientsFromField();

  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > I;
  igl::speye(data.numF, data.numF, I);
  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > QA = data.DDA+lambda*I;
  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > fA = (-2*lambda*Acoeff).sparseView();

  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > QB = data.DDB+lambda*I;
  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > fB = (-2*lambda*I*Bcoeff).sparseView();

  if(doHardConstraints)
  {
    minQuadWithKnownMini(QA, fA, isConstrained, Ak, Acoeff);
    minQuadWithKnownMini(QB, fB, isConstrained, Bk, Bcoeff);
  }
  else
  {
    Eigen::Matrix<int, Eigen::Dynamic, 1>isknown_; isknown_.setZero(data.numF,1);
    Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> xknown_; xknown_.setZero(0,1);
    minQuadWithKnownMini(QA, fA, isknown_, xknown_, Acoeff);
    minQuadWithKnownMini(QB, fB, isknown_, xknown_, Bcoeff);
  }
  setFieldFromCoefficients();

}


template<typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE void igl::AngleBoundFFSolver<DerivedV, DerivedF, DerivedO>::
setFieldFromCoefficients()
{
  for (int i = 0; i <data.numF; ++i)
  {
    //    poly coefficients: 1, 0, -Acoeff, 0, Bcoeff
    //    matlab code from roots (given there are no trailing zeros in the polynomial coefficients)
    Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> polyCoeff(5,1);
    polyCoeff<<1., 0., -Acoeff(i), 0., Bcoeff(i);

    Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> roots;
    polyRoots<std::complex<typename DerivedV::Scalar>>(polyCoeff,roots);

    std::complex<typename DerivedV::Scalar> u = roots[0];
    int maxi = -1;
    float maxd = -1;
    for (int k =1; k<4; ++k)
    {
      float dist = abs(roots[k]+u);
      if (dist>maxd)
      {
        maxd = dist;
        maxi = k;
      }
    }
    std::complex<typename DerivedV::Scalar> v = roots[maxi];
    pvU(i,0) = real(u); pvU(i,1) = imag(u);
    pvV(i,0) = real(v); pvV(i,1) = imag(v);
  }

}

template<typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE void igl::AngleBoundFFSolver<DerivedV, DerivedF, DerivedO>::
minQuadWithKnownMini(const Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > &Q,
                     const Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > &f,
                     const Eigen::VectorXi isConstrained,
                     const Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> &xknown,
                     Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> &x)
{
  int N = Q.rows();

  int nc = xknown.rows();
  Eigen::VectorXi known; known.setZero(nc,1);
  Eigen::VectorXi unknown; unknown.setZero(N-nc,1);

  int indk = 0, indu = 0;
  for (int i = 0; i<N; ++i)
    if (isConstrained[i])
    {
      known[indk] = i;
      indk++;
    }
    else
    {
      unknown[indu] = i;
      indu++;
    }

  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar>> Quu, Quk;

  igl::slice(Q,unknown, unknown, Quu);
  igl::slice(Q,unknown, known, Quk);


  std::vector<typename Eigen::Triplet<std::complex<typename DerivedV::Scalar> > > tripletList;

  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > fu(N-nc,1);

  igl::slice(f,unknown, Eigen::VectorXi::Zero(1,1), fu);

  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > rhs = (Quk*xknown).sparseView()+.5*fu;

  Eigen::SparseLU< Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar>>> solver;
  solver.compute(-Quu);
  if(solver.info()!=Eigen::Success)
  {
    std::cerr<<"Decomposition failed!"<<std::endl;
    return;
  }
  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar>>  b  = solver.solve(rhs);
  if(solver.info()!=Eigen::Success)
  {
    std::cerr<<"Solving failed!"<<std::endl;
    return;
  }

  indk = 0, indu = 0;
  x.setZero(N,1);
  for (int i = 0; i<N; ++i)
    if (isConstrained[i])
      x[i] = xknown[indk++];
    else
      x[i] = b.coeff(indu++,0);

}


template<typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE bool igl::AngleBoundFFSolver<DerivedV, DerivedF, DerivedO>::
solve(const Eigen::VectorXi &isConstrained,
      const Eigen::PlainObjectBase<DerivedO> &initialSolution,
      Eigen::PlainObjectBase<DerivedO> &output,
      typename DerivedV::Scalar *lambdaOut)
{
  int numConstrained = isConstrained.sum();
  // coefficient values
  Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> Ak, Bk;

  pvU.resize(data.numF,2);
  pvV.resize(data.numF,2);
  for (int fi = 0; fi <data.numF; ++fi)
  {
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> &b1 = data.B1.row(fi);
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> &b2 = data.B2.row(fi);
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> &u3 = initialSolution.block(fi,0,1,3);
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> &v3 = initialSolution.block(fi,3,1,3);
    pvU.row(fi)<< u3.dot(b1), u3.dot(b2);
    pvV.row(fi)<< v3.dot(b1), v3.dot(b2);
  }
  setCoefficientsFromField();
  Ak.resize(numConstrained,1);
  Bk.resize(numConstrained,1);
  int ind = 0;
  for (int i = 0; i <data.numF; ++i)
  {
    if(isConstrained[i])
    {
      Ak(ind) = Acoeff[i];
      Bk(ind) = Bcoeff[i];
      ind ++;
    }
  }



  typename DerivedV::Scalar smoothnessValue;
  int oob;

  smoothnessValue = (Acoeff.adjoint()*data.DDA*Acoeff + Bcoeff.adjoint()*data.DDB*Bcoeff).real()[0];
  printf("\n\nInitial smoothness: %.5g\n",smoothnessValue);
  oob = getNumOutOfBounds();
  printf("\n\nInitial out-of-bounds: %d\n",oob);
  printf(" %d %.5g %d\n",-1, smoothnessValue, oob);

  lambda = lambdaInit;
  for (int iter = 0; iter<maxIter; ++iter)
  {
    printf("\n\n--- Iteration %d ---\n",iter);

    localStep();
    globalStep(isConstrained, Ak, Bk);


    smoothnessValue = (Acoeff.adjoint()*data.DDA*Acoeff + Bcoeff.adjoint()*data.DDB*Bcoeff).real()[0];

    printf("Smoothness: %.5g\n",smoothnessValue);

    oob = getNumOutOfBounds();

    bool stoppingCriterion = (oob == 0) ;
    if (stoppingCriterion)
      break;
    lambda = lambda*lambdaMultFactor;
//    printf(" %d %.5g %d\n",iter, smoothnessValue, oob);

  }

  output.setZero(data.numF,6);
  for (int fi=0; fi<data.numF; ++fi)
  {
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> &b1 = data.B1.row(fi);
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> &b2 = data.B2.row(fi);
    output.block(fi,0, 1, 3) = pvU(fi,0)*b1 + pvU(fi,1)*b2;
    output.block(fi,3, 1, 3) = pvV(fi,0)*b1 + pvV(fi,1)*b2;
  }

  if (lambdaOut)
    *lambdaOut = lambda;


  return (oob==0);
}



template <typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE bool igl::angle_bound_frame_fields(const Eigen::PlainObjectBase<DerivedV> &V,
                                            const Eigen::PlainObjectBase<DerivedF> &F,
                                              const typename DerivedV::Scalar &thetaMin,
                                            const Eigen::VectorXi &isConstrained,
                                            const Eigen::PlainObjectBase<DerivedO> &initialSolution,
                                            Eigen::PlainObjectBase<DerivedO> &output,
                                            int maxIter,
                                            const typename DerivedV::Scalar &lambdaInit,
                                            const typename DerivedV::Scalar &lambdaMultFactor,
                                              const bool doHardConstraints)
{
  igl::AngleBoundFFSolverData<DerivedV, DerivedF> csdata(V, F);
  igl::AngleBoundFFSolver<DerivedV, DerivedF, DerivedO> cs(csdata, thetaMin, maxIter, lambdaInit, lambdaMultFactor, doHardConstraints);
  return (cs.solve(isConstrained, initialSolution, output));
}

template <typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE bool igl::angle_bound_frame_fields(const igl::AngleBoundFFSolverData<DerivedV, DerivedF> &csdata,
                                              const typename DerivedV::Scalar &thetaMin,
                                            const Eigen::VectorXi &isConstrained,
                                            const Eigen::PlainObjectBase<DerivedO> &initialSolution,
                                            Eigen::PlainObjectBase<DerivedO> &output,
                                            int maxIter,
                                            const typename DerivedV::Scalar &lambdaInit,
                                            const typename DerivedV::Scalar &lambdaMultFactor,
                                              const bool doHardConstraints,
                                            typename DerivedV::Scalar *lambdaOut)
{
  igl::AngleBoundFFSolver<DerivedV, DerivedF, DerivedO> cs(csdata, thetaMin, maxIter, lambdaInit, lambdaMultFactor, doHardConstraints);
  return (cs.solve(isConstrained, initialSolution, output, lambdaOut));
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
#endif
