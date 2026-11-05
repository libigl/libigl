// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "exact_ray_triangle_hit_compare.h"
#include "exact_ray_triangle_kernel.h"
#include <type_traits>

template <
  typename DerivedS,
  typename DerivedD,
  typename DerivedA>
IGL_INLINE int igl::exact_ray_triangle_hit_compare(
  const Eigen::MatrixBase<DerivedS> & source,
  const Eigen::MatrixBase<DerivedD> & dir,
  const Eigen::MatrixBase<DerivedA> & A1,
  const Eigen::MatrixBase<DerivedA> & B1,
  const Eigen::MatrixBase<DerivedA> & C1,
  const Eigen::MatrixBase<DerivedA> & A2,
  const Eigen::MatrixBase<DerivedA> & B2,
  const Eigen::MatrixBase<DerivedA> & C2)
{
  static_assert(
    std::is_same<typename DerivedA::Scalar, double>::value ||
    std::is_same<typename DerivedA::Scalar, float>::value,
    "exact_ray_triangle_hit_compare requires float or double coordinates");
  namespace k = igl::exact_ray_triangle_kernel;
  const double O[3]  = {(double)source(0),(double)source(1),(double)source(2)};
  const double d[3]  = {(double)dir(0),(double)dir(1),(double)dir(2)};
  const double a1[3] = {(double)A1(0),(double)A1(1),(double)A1(2)};
  const double b1[3] = {(double)B1(0),(double)B1(1),(double)B1(2)};
  const double c1[3] = {(double)C1(0),(double)C1(1),(double)C1(2)};
  const double a2[3] = {(double)A2(0),(double)A2(1),(double)A2(2)};
  const double b2[3] = {(double)B2(0),(double)B2(1),(double)B2(2)};
  const double c2[3] = {(double)C2(0),(double)C2(1),(double)C2(2)};
  return k::compare(O, d, a1, b1, c1, a2, b2, c2);
}

#ifdef IGL_STATIC_LIBRARY
namespace { namespace erthc_types {
  using V3 = Eigen::Matrix<double,3,1,0,3,1>;   // Vector3d
  using R3 = Eigen::Matrix<double,1,3,1,1,3>;   // RowVector3d
}}
#define IGL_ERTHC(S,D,A) \
  template int igl::exact_ray_triangle_hit_compare<S,D,A>( \
    const Eigen::MatrixBase<S>&, const Eigen::MatrixBase<D>&, \
    const Eigen::MatrixBase<A>&, const Eigen::MatrixBase<A>&, const Eigen::MatrixBase<A>&, \
    const Eigen::MatrixBase<A>&, const Eigen::MatrixBase<A>&, const Eigen::MatrixBase<A>&);
IGL_ERTHC(erthc_types::V3, erthc_types::V3, erthc_types::R3)
IGL_ERTHC(erthc_types::R3, erthc_types::R3, erthc_types::R3)
IGL_ERTHC(erthc_types::V3, erthc_types::V3, erthc_types::V3)
#undef IGL_ERTHC
#endif
