// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "exact_ray_triangle_intersect.h"
#include "exact_ray_triangle_kernel.h"
#include <type_traits>

template <
  typename DerivedS,
  typename DerivedD,
  typename DerivedA>
IGL_INLINE bool igl::exact_ray_triangle_intersect(
  const Eigen::MatrixBase<DerivedS> & source,
  const Eigen::MatrixBase<DerivedD> & dir,
  const Eigen::MatrixBase<DerivedA> & A,
  const Eigen::MatrixBase<DerivedA> & B,
  const Eigen::MatrixBase<DerivedA> & C)
{
  static_assert(
    std::is_same<typename DerivedA::Scalar, double>::value ||
    std::is_same<typename DerivedA::Scalar, float>::value,
    "exact_ray_triangle_intersect requires float or double coordinates");
  namespace k = igl::exact_ray_triangle_kernel;
  const double O[3]  = {(double)source(0),(double)source(1),(double)source(2)};
  const double dd[3] = {(double)dir(0),(double)dir(1),(double)dir(2)};
  const double a[3]  = {(double)A(0),(double)A(1),(double)A(2)};
  const double b[3]  = {(double)B(0),(double)B(1),(double)B(2)};
  const double c[3]  = {(double)C(0),(double)C(1),(double)C(2)};
  return k::intersects(O, dd, a, b, c);
}

#ifdef IGL_STATIC_LIBRARY
namespace { namespace erti_types {
  using V3 = Eigen::Matrix<double,3,1,0,3,1>;   // Vector3d
  using R3 = Eigen::Matrix<double,1,3,1,1,3>;   // RowVector3d
}}
#define IGL_ERTI(S,D,A) \
  template bool igl::exact_ray_triangle_intersect<S,D,A>( \
    const Eigen::MatrixBase<S>&, const Eigen::MatrixBase<D>&, \
    const Eigen::MatrixBase<A>&, const Eigen::MatrixBase<A>&, const Eigen::MatrixBase<A>&);
IGL_ERTI(erti_types::V3, erti_types::V3, erti_types::R3)
IGL_ERTI(erti_types::R3, erti_types::R3, erti_types::R3)
IGL_ERTI(erti_types::V3, erti_types::V3, erti_types::V3)
#undef IGL_ERTI
#endif
