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
#include "exactinit.h"

template <
  typename DerivedS,
  typename DerivedD,
  typename DerivedA>
IGL_INLINE int igl::predicates::exact_ray_triangle_hit_compare(
  const Eigen::MatrixBase<DerivedS> & source,
  const Eigen::MatrixBase<DerivedD> & dir,
  const Eigen::MatrixBase<DerivedA> & A1,
  const Eigen::MatrixBase<DerivedA> & B1,
  const Eigen::MatrixBase<DerivedA> & C1,
  const Eigen::MatrixBase<DerivedA> & A2,
  const Eigen::MatrixBase<DerivedA> & B2,
  const Eigen::MatrixBase<DerivedA> & C2)
{
  igl::predicates::exactinit();
  namespace k = igl::predicates::exact_point_triangle_kernel;
  namespace r = igl::predicates::exact_ray_triangle_kernel;
  return r::compare(k::point(source),k::point(dir),k::point(A1),k::point(B1),k::point(C1),
    k::point(A2),k::point(B2),k::point(C2));
}

#ifdef IGL_STATIC_LIBRARY
namespace { namespace erthc_types {
  using V3f = Eigen::Matrix<float,3,1,0,3,1>;
  using R3f = Eigen::Matrix<float,1,3,1,1,3>;
  using V3 = Eigen::Matrix<double,3,1,0,3,1>;   // Vector3d
  using R3 = Eigen::Matrix<double,1,3,1,1,3>;   // RowVector3d
}}
#define IGL_ERTHC(S,D,A) \
  template int igl::predicates::exact_ray_triangle_hit_compare<S,D,A>( \
    const Eigen::MatrixBase<S>&, const Eigen::MatrixBase<D>&, \
    const Eigen::MatrixBase<A>&, const Eigen::MatrixBase<A>&, const Eigen::MatrixBase<A>&, \
    const Eigen::MatrixBase<A>&, const Eigen::MatrixBase<A>&, const Eigen::MatrixBase<A>&);
IGL_ERTHC(erthc_types::V3, erthc_types::V3, erthc_types::R3)
IGL_ERTHC(erthc_types::R3, erthc_types::R3, erthc_types::R3)
IGL_ERTHC(erthc_types::V3, erthc_types::V3, erthc_types::V3)
IGL_ERTHC(erthc_types::R3f, erthc_types::R3f, erthc_types::R3f)
IGL_ERTHC(erthc_types::V3f, erthc_types::V3f, erthc_types::V3f)
#undef IGL_ERTHC
#endif
