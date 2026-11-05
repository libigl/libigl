// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_EXACT_RAY_TRIANGLE_INTERSECT_H
#define IGL_EXACT_RAY_TRIANGLE_INTERSECT_H
#include "igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  /// Exact 3-way classification of a ray against a triangle, returned by
  /// exact_ray_triangle_intersect.
  enum class RayTriangleIntersection : int
  {
    /// The ray does not intersect the triangle.
    None = 0,
    /// The ray transversally intersects the triangle (a clean point hit, t ≥ 0).
    Hit = 1,
    /// The ray lies exactly in the triangle's supporting plane. The intersection,
    /// if any, is a segment rather than a point; this routine does not resolve
    /// whether/where the coplanar ray overlaps the triangle.
    Coplanar = 2
  };

  /// Exactly classify the ray {source + t·dir : t ≥ 0} against the triangle
  /// (A,B,C) using only exact sign predicates (adaptive-precision expansion
  /// arithmetic). No intersection point is constructed and no rounded
  /// intermediate quantity (such as source+dir) is ever formed: the sign tests
  /// refer to the exact mathematical ray defined by the floating point inputs
  /// `source` and `dir`. There are no tolerances or fudge factors — every branch
  /// (including the coplanar case) is decided by an exact determinant sign.
  ///
  /// Transverse intersections (including grazing an edge or vertex) return `Hit`.
  /// If the ray direction is exactly parallel to the triangle's plane and the
  /// origin lies exactly in that plane, the configuration is `Coplanar`;
  /// otherwise (parallel but offset) it is `None`.
  ///
  /// The exact arithmetic assumes IEEE-754 double; coordinates must be `float` or
  /// `double` (both promote to double exactly, so the classification is exact
  /// with respect to the given coordinates).
  ///
  /// @param[in] source  3-vector origin of ray
  /// @param[in] dir     3-vector direction of ray
  /// @param[in] A       3-vector first triangle corner
  /// @param[in] B       3-vector second triangle corner
  /// @param[in] C       3-vector third triangle corner
  /// @return None, Hit, or Coplanar (see RayTriangleIntersection)
  ///
  /// \see exact_ray_mesh_intersect, exact_ray_triangle_hit_compare
  template <
    typename DerivedS,
    typename DerivedD,
    typename DerivedA>
  IGL_INLINE RayTriangleIntersection exact_ray_triangle_intersect(
    const Eigen::MatrixBase<DerivedS> & source,
    const Eigen::MatrixBase<DerivedD> & dir,
    const Eigen::MatrixBase<DerivedA> & A,
    const Eigen::MatrixBase<DerivedA> & B,
    const Eigen::MatrixBase<DerivedA> & C);
}

#ifndef IGL_STATIC_LIBRARY
#  include "exact_ray_triangle_intersect.cpp"
#endif
#endif
