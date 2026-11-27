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
#include "exactinit.h"

template <
  typename DerivedS,
  typename DerivedD,
  typename DerivedA>
IGL_INLINE igl::predicates::RayTriangleIntersection igl::predicates::exact_ray_triangle_intersect(
  const Eigen::MatrixBase<DerivedS> & source,
  const Eigen::MatrixBase<DerivedD> & dir,
  const Eigen::MatrixBase<DerivedA> & A,
  const Eigen::MatrixBase<DerivedA> & B,
  const Eigen::MatrixBase<DerivedA> & C)
{
  igl::predicates::exactinit();
  namespace k = igl::predicates::exact_point_triangle_kernel;
  namespace r = igl::predicates::exact_ray_triangle_kernel;
  return static_cast<igl::predicates::RayTriangleIntersection>(r::classify(k::point(source),k::point(dir),k::point(A),k::point(B),k::point(C)));
}

#ifdef IGL_STATIC_LIBRARY
namespace { namespace erti_types {
  using V3f = Eigen::Matrix<float,3,1,0,3,1>;
  using R3f = Eigen::Matrix<float,1,3,1,1,3>;
  using V3 = Eigen::Matrix<double,3,1,0,3,1>;   // Vector3d
  using R3 = Eigen::Matrix<double,1,3,1,1,3>;   // RowVector3d
}}
#define IGL_ERTI(S,D,A) \
  template igl::predicates::RayTriangleIntersection igl::predicates::exact_ray_triangle_intersect<S,D,A>( \
    const Eigen::MatrixBase<S>&, const Eigen::MatrixBase<D>&, \
    const Eigen::MatrixBase<A>&, const Eigen::MatrixBase<A>&, const Eigen::MatrixBase<A>&);
IGL_ERTI(erti_types::V3, erti_types::V3, erti_types::R3)
IGL_ERTI(erti_types::R3, erti_types::R3, erti_types::R3)
IGL_ERTI(erti_types::V3, erti_types::V3, erti_types::V3)
IGL_ERTI(erti_types::R3f, erti_types::R3f, erti_types::R3f)
IGL_ERTI(erti_types::V3f, erti_types::V3f, erti_types::V3f)
#undef IGL_ERTI
#endif
