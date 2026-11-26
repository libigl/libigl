// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "point_mesh_closest_point.h"
#include "point_mesh_closest_point_kernel.h"

template<typename DerivedQ,typename DerivedV,typename DerivedF,
  typename DerivedB,typename DerivedL,typename DerivedBary>
IGL_INLINE double igl::predicates::point_mesh_closest_point(
  const Eigen::MatrixBase<DerivedQ> &query,
  const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
  const Eigen::MatrixBase<DerivedL> &leaf,
  Eigen::Index &I,Eigen::PlainObjectBase<DerivedBary> &bary)
{
  const auto closest=point_mesh_closest_point_kernel::closest_point(query,V,F,B1,B2,leaf);
  I=closest.triangle;
  point_mesh_closest_point_kernel::barycentrics(closest,bary);
  return closest.distance;
}

#ifdef IGL_STATIC_LIBRARY
#define IGL_POINT_MESH_CLOSEST_POINT(Q) \
template double igl::predicates::point_mesh_closest_point< \
  Q,Eigen::MatrixXd,Eigen::MatrixXi,Eigen::MatrixXd,Eigen::VectorXi,Eigen::RowVector3d>( \
  const Eigen::MatrixBase<Q>&,const Eigen::MatrixBase<Eigen::MatrixXd>&, \
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXd>&, \
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::VectorXi>&, \
  Eigen::Index&,Eigen::PlainObjectBase<Eigen::RowVector3d>&);
IGL_POINT_MESH_CLOSEST_POINT(Eigen::RowVector3d)
IGL_POINT_MESH_CLOSEST_POINT(Eigen::Vector3d)
#undef IGL_POINT_MESH_CLOSEST_POINT
#endif
