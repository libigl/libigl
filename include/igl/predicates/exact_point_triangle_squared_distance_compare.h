// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PREDICATES_EXACT_POINT_TRIANGLE_SQUARED_DISTANCE_COMPARE_H
#define IGL_PREDICATES_EXACT_POINT_TRIANGLE_SQUARED_DISTANCE_COMPARE_H
#include "../igl_inline.h"
#include <Eigen/Core>
namespace igl { namespace predicates {
  /// Compare distances from a query to two closed triangles, exactly.
  /// Closest points and squared distances are not constructed. Degenerate
  /// triangles are supported with the same contract as closest_feature.
  /// Inputs must be finite float or double coordinates. Link igl::predicates
  /// with its double backend. Correctness assumes round-to-nearest without
  /// reassociation, contraction, or intermediate overflow/underflow, as for
  /// Shewchuk's predicates. No expansion survives the call.
  /// @param[in] query  Query 3-vector
  /// @param[in] a1,b1,c1,a2,b2,c2  First and second triangle corners
  /// @return -1 if the first is closer, 0 for a tie, +1 if farther
  template<typename DerivedQ,typename DerivedA>
  IGL_INLINE int exact_point_triangle_squared_distance_compare(
    const Eigen::MatrixBase<DerivedQ> &query,
    const Eigen::MatrixBase<DerivedA> &a1,
    const Eigen::MatrixBase<DerivedA> &b1,
    const Eigen::MatrixBase<DerivedA> &c1,
    const Eigen::MatrixBase<DerivedA> &a2,
    const Eigen::MatrixBase<DerivedA> &b2,
    const Eigen::MatrixBase<DerivedA> &c2);
}}
#ifndef IGL_STATIC_LIBRARY
#  include "exact_point_triangle_squared_distance_compare.cpp"
#endif
#endif
