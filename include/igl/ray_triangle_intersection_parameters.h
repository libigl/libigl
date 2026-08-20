// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_RAY_TRIANGLE_INTERSECTION_PARAMETERS_H
#define IGL_RAY_TRIANGLE_INTERSECTION_PARAMETERS_H
#include "igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  /// Solve for the parameters (t,u,v) where the supporting line of the ray
  /// {source + t·dir} meets the plane of triangle (A,B,C), via Möller–Trumbore:
  ///
  ///   source + t·dir = A + u·(B-A) + v·(C-A)
  ///
  /// This is an ordinary (inexact) floating point solve. Unlike
  /// igl::ray_triangle_intersect it performs NO culling: it does not check the
  /// barycentric bounds (0 ≤ u, v, u+v ≤ 1) or the sign of t. It is intended for
  /// constructing an approximate hit record after intersection has already been
  /// classified (e.g. by exact_ray_triangle_intersect), or as a building block.
  ///
  /// @param[in] source  3-vector origin of ray
  /// @param[in] dir     3-vector direction of ray
  /// @param[in] A       3-vector first triangle corner
  /// @param[in] B       3-vector second triangle corner
  /// @param[in] C       3-vector third triangle corner
  /// @param[out] t  ray parameter of the plane intersection
  /// @param[out] u  barycentric coordinate along edge (B-A)
  /// @param[out] v  barycentric coordinate along edge (C-A)
  /// @return false if the ray is parallel to the triangle's plane (t,u,v set to 0)
  ///
  /// \see ray_triangle_intersect, exact_ray_triangle_intersect
  template <
    typename DerivedS,
    typename DerivedD,
    typename DerivedA,
    typename Scalar>
  IGL_INLINE bool ray_triangle_intersection_parameters(
    const Eigen::MatrixBase<DerivedS> & source,
    const Eigen::MatrixBase<DerivedD> & dir,
    const Eigen::MatrixBase<DerivedA> & A,
    const Eigen::MatrixBase<DerivedA> & B,
    const Eigen::MatrixBase<DerivedA> & C,
    Scalar & t,
    Scalar & u,
    Scalar & v);
}

#ifndef IGL_STATIC_LIBRARY
#  include "ray_triangle_intersection_parameters.cpp"
#endif
#endif
