// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PREDICATES_EXACT_RAY_TRIANGLE_INTERSECT_H
#define IGL_PREDICATES_EXACT_RAY_TRIANGLE_INTERSECT_H
#include "../igl_inline.h"
#include <Eigen/Core>

namespace igl
{
namespace predicates
{
  /// Exact 3-way classification of a ray against a triangle, returned by
  /// exact_ray_triangle_intersect.
  enum class RayTriangleIntersection : int
  {
    /// The ray does not intersect the triangle.
    None = 0,
    /// The ray transversally intersects the triangle (a clean point hit, t ≥ 0).
    Hit = 1,
    /// The ray lies exactly in the triangle's supporting plane. This routine
    /// does not resolve whether the coplanar intersection is empty, a point,
    /// or a segment.
    Coplanar = 2
  };

  /// Exactly classify the ray {source + t·dir : t ≥ 0} against the triangle
  /// (A,B,C) using only exact sign predicates (adaptive-precision expansion
  /// arithmetic). No intersection point or rounded second ray point
  /// (such as source+dir) is constructed: the sign tests refer to the exact
  /// mathematical ray defined by the floating point inputs
  /// `source` and `dir`. There are no tolerances or fudge factors — every branch
  /// (including the coplanar case) is decided by an exact determinant sign.
  ///
  /// Transverse intersections (including grazing an edge or vertex) return `Hit`.
  /// If the ray direction is exactly parallel to the triangle's plane and the
  /// origin lies exactly in that plane, the configuration is `Coplanar`;
  /// otherwise (parallel but offset) it is `None`.
  ///
  /// The exact arithmetic assumes IEEE-754 double; coordinates must be finite
  /// `float` or `double` (both promote to double exactly, so classification is exact
  /// with respect to the given coordinates).
  ///
  /// Link against igl::predicates with its double-precision backend. Correctness
  /// assumes round-to-nearest with gradual underflow, without reassociation or
  /// contraction. Exact refinement rescales positions and the direction by
  /// separate positive powers of two. It throws std::range_error if an exponent
  /// span cannot fit safely. A sufficient range for both ray predicates is:
  /// after rescaling, every nonzero component's lowest set bit is at least
  /// 2^-166 and its magnitude is less than 2^151. Float inputs promote exactly
  /// to double. No expansion is retained between calls.
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
}

#ifndef IGL_STATIC_LIBRARY
#  include "exact_ray_triangle_intersect.cpp"
#endif
#endif
