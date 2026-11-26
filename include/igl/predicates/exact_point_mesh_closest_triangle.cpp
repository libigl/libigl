// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "exact_point_mesh_closest_triangle.h"
#include "exact_point_mesh_closest_triangle_kernel.h"

template<typename DerivedQ,typename DerivedV,typename DerivedF,
  typename DerivedB,typename DerivedL>
IGL_INLINE Eigen::Index igl::predicates::exact_point_mesh_closest_triangle(
  const Eigen::MatrixBase<DerivedQ> &query,
  const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
  const Eigen::MatrixBase<DerivedL> &leaf)
{
  return exact_point_mesh_closest_triangle_kernel::closest_feature(query,V,F,B1,B2,leaf).triangle;
}

#ifdef IGL_STATIC_LIBRARY
template Eigen::Index igl::predicates::exact_point_mesh_closest_triangle<
  Eigen::RowVector3d,Eigen::MatrixXd,Eigen::MatrixXi,Eigen::MatrixXd,Eigen::VectorXi>(
  const Eigen::MatrixBase<Eigen::RowVector3d>&,const Eigen::MatrixBase<Eigen::MatrixXd>&,
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXd>&,
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::VectorXi>&);
template Eigen::Index igl::predicates::exact_point_mesh_closest_triangle<
  Eigen::Vector3d,Eigen::MatrixXd,Eigen::MatrixXi,Eigen::MatrixXd,Eigen::VectorXi>(
  const Eigen::MatrixBase<Eigen::Vector3d>&,const Eigen::MatrixBase<Eigen::MatrixXd>&,
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXd>&,
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::VectorXi>&);
#endif
