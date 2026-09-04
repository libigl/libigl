// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_EXACT_RAY_TRIANGLE_HIT_COMPARE_H
#define IGL_EXACT_RAY_TRIANGLE_HIT_COMPARE_H
#include "igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  /// Exactly compare the (unconstructed) hit distances of two triangles along a
  /// common ray. Both triangles are assumed to be transversally hit by the ray
  /// (see exact_ray_triangle_intersect). The comparison is evaluated as the exact
  /// sign of a difference of products of determinants computed directly from the
  /// original floating point coordinates — no ratio t = -a/b is ever formed, so no
  /// exact construction is exposed.
  ///
  /// Coordinates must be `float` or `double` (both promote to double exactly).
  ///
  /// @param[in] source  3-vector origin of ray
  /// @param[in] dir     3-vector direction of ray
  /// @param[in] A1,B1,C1  corners of first triangle
  /// @param[in] A2,B2,C2  corners of second triangle
  /// @return -1 if t1 < t2, 0 if t1 == t2, +1 if t1 > t2
  ///
  /// \see exact_ray_triangle_intersect, exact_ray_mesh_intersect
  template <
    typename DerivedS,
    typename DerivedD,
    typename DerivedA>
  IGL_INLINE int exact_ray_triangle_hit_compare(
    const Eigen::MatrixBase<DerivedS> & source,
    const Eigen::MatrixBase<DerivedD> & dir,
    const Eigen::MatrixBase<DerivedA> & A1,
    const Eigen::MatrixBase<DerivedA> & B1,
    const Eigen::MatrixBase<DerivedA> & C1,
    const Eigen::MatrixBase<DerivedA> & A2,
    const Eigen::MatrixBase<DerivedA> & B2,
    const Eigen::MatrixBase<DerivedA> & C2);
}

#ifndef IGL_STATIC_LIBRARY
#  include "exact_ray_triangle_hit_compare.cpp"
#endif
#endif
