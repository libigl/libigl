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
  /// Exactly determine whether the ray {source + t·dir : t ≥ 0} intersects the
  /// triangle (A,B,C) using only exact sign predicates (adaptive-precision
  /// expansion arithmetic). No intersection point is constructed and no rounded
  /// intermediate quantity (such as source+dir) is ever formed: the sign tests
  /// refer to the exact mathematical ray defined by the floating point inputs
  /// `source` and `dir`.
  ///
  /// Transverse intersections (including grazing an edge or vertex) are decided
  /// exactly. A ray whose direction is exactly parallel to the triangle's
  /// supporting plane (the denominator determinant det(dir,B-A,C-A) is exactly
  /// zero) is a measure-zero degeneracy whose intersection is generally a
  /// t-interval rather than a point; this routine reports such coplanar
  /// configurations as a non-hit.
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
  /// @return true if the ray transversally intersects the triangle with t ≥ 0
  ///
  /// \see exact_ray_mesh_intersect, exact_ray_triangle_hit_compare
  template <
    typename DerivedS,
    typename DerivedD,
    typename DerivedA>
  IGL_INLINE bool exact_ray_triangle_intersect(
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
