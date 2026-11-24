// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "exact_point_triangle_squared_distance_compare.h"
#include "exact_point_triangle_kernel.h"
#include "exactinit.h"
template<typename DerivedQ,typename DerivedA>
IGL_INLINE int igl::predicates::exact_point_triangle_squared_distance_compare(
  const Eigen::MatrixBase<DerivedQ> &query,
  const Eigen::MatrixBase<DerivedA> &a1,
  const Eigen::MatrixBase<DerivedA> &b1,
  const Eigen::MatrixBase<DerivedA> &c1,
  const Eigen::MatrixBase<DerivedA> &a2,
  const Eigen::MatrixBase<DerivedA> &b2,
  const Eigen::MatrixBase<DerivedA> &c2)
{
  igl::predicates::exactinit();
  namespace k=exact_point_triangle_kernel;
  return k::compare(k::point(query),{{k::point(a1),k::point(b1),k::point(c1)}},{{k::point(a2),k::point(b2),k::point(c2)}});
}

#ifdef IGL_STATIC_LIBRARY
template int igl::predicates::exact_point_triangle_squared_distance_compare<Eigen::Matrix<double,1,3>,Eigen::Matrix<double,1,3>>(const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &);
template int igl::predicates::exact_point_triangle_squared_distance_compare<Eigen::Matrix<double,3,1>,Eigen::Matrix<double,3,1>>(const Eigen::MatrixBase<Eigen::Matrix<double,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<double,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<double,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<double,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<double,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<double,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<double,3,1>> &);
template int igl::predicates::exact_point_triangle_squared_distance_compare<Eigen::Matrix<double,3,1>,Eigen::Matrix<double,1,3>>(const Eigen::MatrixBase<Eigen::Matrix<double,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &);
template int igl::predicates::exact_point_triangle_squared_distance_compare<Eigen::Matrix<float,1,3>,Eigen::Matrix<float,1,3>>(const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &);
template int igl::predicates::exact_point_triangle_squared_distance_compare<Eigen::Matrix<float,3,1>,Eigen::Matrix<float,3,1>>(const Eigen::MatrixBase<Eigen::Matrix<float,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<float,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<float,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<float,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<float,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<float,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<float,3,1>> &);
template int igl::predicates::exact_point_triangle_squared_distance_compare<Eigen::Matrix<float,3,1>,Eigen::Matrix<float,1,3>>(const Eigen::MatrixBase<Eigen::Matrix<float,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &);
#endif
