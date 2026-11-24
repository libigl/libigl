// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "exact_point_triangle_closest_feature.h"
#include "exact_point_triangle_kernel.h"
#include "exactinit.h"
template<typename DerivedQ,typename DerivedA>
IGL_INLINE igl::predicates::ClosestPointFeature igl::predicates::exact_point_triangle_closest_feature(
  const Eigen::MatrixBase<DerivedQ> &query,
  const Eigen::MatrixBase<DerivedA> &a,
  const Eigen::MatrixBase<DerivedA> &b,
  const Eigen::MatrixBase<DerivedA> &c,int &feature_id)
{
  igl::predicates::exactinit();
  namespace k=exact_point_triangle_kernel;
  return static_cast<ClosestPointFeature>(k::feature(k::point(query),{{k::point(a),k::point(b),k::point(c)}},feature_id));
}

#ifdef IGL_STATIC_LIBRARY
template igl::predicates::ClosestPointFeature igl::predicates::exact_point_triangle_closest_feature<Eigen::Matrix<double,1,3>,Eigen::Matrix<double,1,3>>(const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, int &);
template igl::predicates::ClosestPointFeature igl::predicates::exact_point_triangle_closest_feature<Eigen::Matrix<double,3,1>,Eigen::Matrix<double,3,1>>(const Eigen::MatrixBase<Eigen::Matrix<double,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<double,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<double,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<double,3,1>> &, int &);
template igl::predicates::ClosestPointFeature igl::predicates::exact_point_triangle_closest_feature<Eigen::Matrix<double,3,1>,Eigen::Matrix<double,1,3>>(const Eigen::MatrixBase<Eigen::Matrix<double,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<double,1,3>> &, int &);
template igl::predicates::ClosestPointFeature igl::predicates::exact_point_triangle_closest_feature<Eigen::Matrix<float,1,3>,Eigen::Matrix<float,1,3>>(const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, int &);
template igl::predicates::ClosestPointFeature igl::predicates::exact_point_triangle_closest_feature<Eigen::Matrix<float,3,1>,Eigen::Matrix<float,3,1>>(const Eigen::MatrixBase<Eigen::Matrix<float,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<float,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<float,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<float,3,1>> &, int &);
template igl::predicates::ClosestPointFeature igl::predicates::exact_point_triangle_closest_feature<Eigen::Matrix<float,3,1>,Eigen::Matrix<float,1,3>>(const Eigen::MatrixBase<Eigen::Matrix<float,3,1>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, const Eigen::MatrixBase<Eigen::Matrix<float,1,3>> &, int &);
#endif
