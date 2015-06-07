// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "signed_distance.h"
#include "per_vertex_normals.h"
#include "per_edge_normals.h"
#include "per_face_normals.h"
#include "get_seconds.h"
#include "point_mesh_squared_distance.h"


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
  assert(V.cols() == 3 && "V should have 3d positions");
  assert(P.cols() == 3 && "P should have 3d positions");
  // Only unsigned distance is supported for non-triangles
  if(sign_type != SIGNED_DISTANCE_TYPE_UNSIGNED)
  {
    assert(F.cols() == 3 && "F should have triangles");
  }

  // Prepare distance computation
  AABB<MatrixXd,3> tree;
  tree.init(V,F);

  Eigen::MatrixXd FN,VN,EN;
  Eigen::MatrixXi E;
  Eigen::VectorXi EMAP;
  WindingNumberAABB<Eigen::Vector3d> hier;
  switch(sign_type)
  {
    default:
      assert(false && "Unknown SignedDistanceType");
    case SIGNED_DISTANCE_TYPE_UNSIGNED:
      // do nothing
      break;
    case SIGNED_DISTANCE_TYPE_DEFAULT:
    case SIGNED_DISTANCE_TYPE_WINDING_NUMBER:
      hier.set_mesh(V,F);
      hier.grow();
      break;
    case SIGNED_DISTANCE_TYPE_PSEUDONORMAL:
      // "Signed Distance Computation Using the Angle Weighted Pseudonormal"
      // [Bærentzen & Aanæs 2005]
      per_face_normals(V,F,FN);
      per_vertex_normals(V,F,PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE,FN,VN);
      per_edge_normals(
        V,F,PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM,FN,EN,E,EMAP);
      N.resize(P.rows(),3);
      break;
  }

  S.resize(P.rows(),1);
  I.resize(P.rows(),1);
  C.resize(P.rows(),3);
  for(int p = 0;p<P.rows();p++)
  {
    const RowVector3d q = P.row(p);
    double s,sqrd;
    RowVector3d c;
    int i=-1;
    switch(sign_type)
    {
      default:
        assert(false && "Unknown SignedDistanceType");
      case SIGNED_DISTANCE_TYPE_UNSIGNED:
        s = 1.;
        sqrd = tree.squared_distance(V,F,q,i,c);
        break;
      case SIGNED_DISTANCE_TYPE_DEFAULT:
      case SIGNED_DISTANCE_TYPE_WINDING_NUMBER:
        signed_distance_winding_number(tree,V,F,hier,q,s,sqrd,i,c);
        break;
      case SIGNED_DISTANCE_TYPE_PSEUDONORMAL:
      {
        Eigen::RowVector3d n;
        signed_distance_pseudonormal(tree,V,F,FN,VN,EN,EMAP,q,s,sqrd,i,c,n);
        N.row(p) = n;
        break;
      }
    }
    I(p) = i;
    S(p) = s*sqrt(sqrd);
    C.row(p) = c;
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
  const auto & qc = q-c;
  RowVector3d b;
  AABB<Eigen::MatrixXd,3>::barycentric_coordinates(
    c,V.row(F(f,0)),V.row(F(f,1)),V.row(F(f,2)),b);
  // Determine which normal to use
  const double epsilon = 1e-12;
  const int type = (b.array()<=epsilon).cast<int>().sum();
  switch(type)
  {
    case 2:
      // Find vertex
      for(int x = 0;x<3;x++)
      {
        if(b(x)>epsilon)
        {
          n = VN.row(F(f,x));
          break;
        }
      }
      break;
    case 1:
      // Find edge
      for(int x = 0;x<3;x++)
      {
        if(b(x)<=epsilon)
        {
          n = EN.row(EMAP(F.rows()*x+f));
          break;
        }
      }
      break;
    default:
      assert(false && "all barycentric coords zero.");
    case 0:
      n = FN.row(f);
      break;
  }
  s = (qc.dot(n) >= 0 ? 1. : -1.);
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
  const igl::WindingNumberAABB<Eigen::Vector3d> & hier,
  const Eigen::RowVector3d & q,
  double & s,
  double & sqrd,
  int & i,
  Eigen::RowVector3d & c)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  sqrd = tree.squared_distance(V,F,q,i,c);
  const double w = hier.winding_number(q.transpose());
  s = 1.-2.*w;
}

#ifdef IGL_STATIC_LIBRARY
// This template is necessary for the others to compile with clang
// http://stackoverflow.com/questions/27748442/is-clangs-c11-support-reliable
#endif
