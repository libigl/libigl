// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "signed_distance.h"
#include "get_seconds.h"
#include "per_edge_normals.h"
#include "per_face_normals.h"
#include "per_vertex_normals.h"
#include "point_mesh_squared_distance.h"
#include "pseudonormal_test.h"


IGL_INLINE void igl::signed_distance(
  const Eigen::MatrixXd & P,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const SignedDistanceType sign_type,
  Eigen::VectorXd & S,
  Eigen::VectorXi & I,
  Eigen::MatrixXd & C,
  Eigen::MatrixXd & N)
{
  using namespace Eigen;
  using namespace std;
  const int dim = V.cols();
  assert((V.cols() == 3||V.cols() == 2) && "V should have 3d or 2d positions");
  assert((P.cols() == 3||P.cols() == 2) && "P should have 3d or 2d positions");
  assert(V.cols() == P.cols() && "V should have same dimension as P");
  // Only unsigned distance is supported for non-triangles
  if(sign_type != SIGNED_DISTANCE_TYPE_UNSIGNED)
  {
    assert(F.cols() == dim && "F should have co-dimension 0 simplices");
  }

  // Prepare distance computation
  AABB<MatrixXd,3> tree3;
  AABB<MatrixXd,2> tree2;
  switch(dim)
  {
    default:
    case 3:
      tree3.init(V,F);
      break;
    case 2:
      tree2.init(V,F);
      break;
  }

  Eigen::MatrixXd FN,VN,EN;
  Eigen::MatrixXi E;
  Eigen::VectorXi EMAP;
  WindingNumberAABB<Eigen::Vector3d> hier3;
  switch(sign_type)
  {
    default:
      assert(false && "Unknown SignedDistanceType");
    case SIGNED_DISTANCE_TYPE_UNSIGNED:
      // do nothing
      break;
    case SIGNED_DISTANCE_TYPE_DEFAULT:
    case SIGNED_DISTANCE_TYPE_WINDING_NUMBER:
      switch(dim)
      {
        default:
        case 3:
          hier3.set_mesh(V,F);
          hier3.grow();
          break;
        case 2:
          // no precomp, no hierarchy
          break;
      }
      break;
    case SIGNED_DISTANCE_TYPE_PSEUDONORMAL:
      switch(dim)
      {
        default:
        case 3:
          // "Signed Distance Computation Using the Angle Weighted Pseudonormal"
          // [Bærentzen & Aanæs 2005]
          per_face_normals(V,F,FN);
          per_vertex_normals(V,F,PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE,FN,VN);
          per_edge_normals(
            V,F,PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM,FN,EN,E,EMAP);
          break;
        case 2:
          FN.resize(F.rows(),2);
          VN = MatrixXd::Zero(V.rows(),2);
          for(int e = 0;e<F.rows();e++)
          {
            // rotate edge vector
            FN(e,0) = -(V(F(e,1),1)-V(F(e,0),1));
            FN(e,1) =  (V(F(e,1),0)-V(F(e,0),0));
            FN.row(e).normalize();
            // add to vertex normal
            VN.row(F(e,1)) += FN.row(e);
            VN.row(F(e,0)) += FN.row(e);
          }
          // normalize to average
          VN.rowwise().normalize();
          break;
      }
      N.resize(P.rows(),dim);
      break;
  }

  S.resize(P.rows(),1);
  I.resize(P.rows(),1);
  C.resize(P.rows(),dim);
  for(int p = 0;p<P.rows();p++)
  {
    RowVector3d q3;
    RowVector2d q2;
    switch(P.cols())
    {
      default:
      case 3:
        q3 = P.row(p);
        break;
      case 2:
        q2 = P.row(p);
        break;
    }
    double s,sqrd;
    RowVectorXd c;
    RowVector3d c3;
    RowVector2d c2;
    int i=-1;
    switch(sign_type)
    {
      default:
        assert(false && "Unknown SignedDistanceType");
      case SIGNED_DISTANCE_TYPE_UNSIGNED:
        s = 1.;
        sqrd = dim==3?
          tree3.squared_distance(V,F,q3,i,c3):
          tree2.squared_distance(V,F,q2,i,c2);
        break;
      case SIGNED_DISTANCE_TYPE_DEFAULT:
      case SIGNED_DISTANCE_TYPE_WINDING_NUMBER:
        dim==3 ? 
          signed_distance_winding_number(tree3,V,F,hier3,q3,s,sqrd,i,c3):
          signed_distance_winding_number(tree2,V,F,q2,s,sqrd,i,c2);
        break;
      case SIGNED_DISTANCE_TYPE_PSEUDONORMAL:
      {
        RowVector3d n3;
        RowVector2d n2;
        dim==3 ?
          signed_distance_pseudonormal(tree3,V,F,FN,VN,EN,EMAP,q3,s,sqrd,i,c3,n3):
          signed_distance_pseudonormal(tree2,V,F,FN,VN,q2,s,sqrd,i,c2,n2);
        Eigen::RowVectorXd n;
        (dim==3 ? n = n3 : n = n2);
        N.row(p) = n;
        break;
      }
    }
    I(p) = i;
    S(p) = s*sqrt(sqrd);
    C.row(p) = (dim==3 ? c=c3 : c=c2);
  }
}


IGL_INLINE double igl::signed_distance_pseudonormal(
  const AABB<Eigen::MatrixXd,3> & tree,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & FN,
  const Eigen::MatrixXd & VN,
  const Eigen::MatrixXd & EN,
  const Eigen::VectorXi & EMAP,
  const Eigen::RowVector3d & q)
{
  double s,sqrd;
  Eigen::RowVector3d n,c;
  int i = -1;
  signed_distance_pseudonormal(tree,V,F,FN,VN,EN,EMAP,q,s,sqrd,i,c,n);
  return s*sqrt(sqrd);
}

IGL_INLINE void igl::signed_distance_pseudonormal(
  const Eigen::MatrixXd & P,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const AABB<Eigen::MatrixXd,3> & tree,
  const Eigen::MatrixXd & FN,
  const Eigen::MatrixXd & VN,
  const Eigen::MatrixXd & EN,
  const Eigen::VectorXi & EMAP,
  Eigen::VectorXd & S,
  Eigen::VectorXi & I,
  Eigen::MatrixXd & C,
  Eigen::MatrixXd & N)
{
  using namespace Eigen;
  const size_t np = P.rows();
  S.resize(np,1);
  I.resize(np,1);
  N.resize(np,3);
  C.resize(np,3);
# pragma omp parallel for if(np>1000)
  for(size_t p = 0;p<np;p++)
  {
    double s,sqrd;
    RowVector3d n,c;
    int i = -1;
    RowVector3d q = P.row(p);
    signed_distance_pseudonormal(tree,V,F,FN,VN,EN,EMAP,q,s,sqrd,i,c,n);
    S(p) = s*sqrt(sqrd);
    I(p) = i;
    N.row(p) = n;
    C.row(p) = c;
  }
//  igl::AABB<MatrixXd,3> tree_P;
//  MatrixXi J = VectorXi::LinSpaced(P.rows(),0,P.rows()-1);
//  tree_P.init(P,J);
//  tree.squared_distance(V,F,tree_P,P,J,S,I,C);
//# pragma omp parallel for if(np>1000)
//  for(size_t p = 0;p<np;p++)
//  {
//    RowVector3d c = C.row(p);
//    RowVector3d q = P.row(p);
//    const int f = I(p);
//    double s;
//    RowVector3d n;
//    pseudonormal_test(V,F,FN,VN,EN,EMAP,q,f,c,s,n);
//    N.row(p) = n;
//    S(p) = s*sqrt(S(p));
//  }

}

IGL_INLINE void igl::signed_distance_pseudonormal(
  const AABB<Eigen::MatrixXd,3> & tree,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & FN,
  const Eigen::MatrixXd & VN,
  const Eigen::MatrixXd & EN,
  const Eigen::VectorXi & EMAP,
  const Eigen::RowVector3d & q,
  double & s,
  double & sqrd,
  int & f,
  Eigen::RowVector3d & c,
  Eigen::RowVector3d & n)
{
  using namespace Eigen;
  using namespace std;
  sqrd = tree.squared_distance(V,F,q,f,c);
  pseudonormal_test(V,F,FN,VN,EN,EMAP,q,f,c,s,n);
}

IGL_INLINE void igl::signed_distance_pseudonormal(
  const AABB<Eigen::MatrixXd,2> & tree,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & FN,
  const Eigen::MatrixXd & VN,
  const Eigen::RowVector2d & q,
  double & s,
  double & sqrd,
  int & f,
  Eigen::RowVector2d & c,
  Eigen::RowVector2d & n)
{
  using namespace Eigen;
  using namespace std;
  sqrd = tree.squared_distance(V,F,q,f,c);
  pseudonormal_test(V,F,FN,VN,q,f,c,s,n);
}

IGL_INLINE double igl::signed_distance_winding_number(
  const AABB<Eigen::MatrixXd,3> & tree,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const igl::WindingNumberAABB<Eigen::Vector3d> & hier,
  const Eigen::RowVector3d & q)
{
  double s,sqrd;
  Eigen::RowVector3d c;
  int i=-1;
  signed_distance_winding_number(tree,V,F,hier,q,s,sqrd,i,c);
  return s*sqrt(sqrd);
}


IGL_INLINE void igl::signed_distance_winding_number(
  const AABB<Eigen::MatrixXd,3> & tree,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const igl::WindingNumberAABB<Eigen::Matrix<double,3,1> > & hier,
  const Eigen::Matrix<double,1,3> & q,
  double & s,
  double & sqrd,
  int & i,
  Eigen::Matrix<double,1,3> & c)
{
  using namespace Eigen;
  using namespace std;
  sqrd = tree.squared_distance(V,F,q,i,c);
  const double w = hier.winding_number(q.transpose());
  s = 1.-2.*w;
}

IGL_INLINE void igl::signed_distance_winding_number(
  const AABB<Eigen::MatrixXd,2> & tree,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::Matrix<double,1,2> & q,
  double & s,
  double & sqrd,
  int & i,
  Eigen::Matrix<double,1,2> & c)
{
  using namespace Eigen;
  using namespace std;
  sqrd = tree.squared_distance(V,F,q,i,c);
  double w;
  winding_number_2(V.data(), V.rows(), F.data(), F.rows(), q.data(), 1, &w);
  s = 1.-2.*w;

}

#ifdef IGL_STATIC_LIBRARY
// This template is necessary for the others to compile with clang
// http://stackoverflow.com/questions/27748442/is-clangs-c11-support-reliable
#endif
