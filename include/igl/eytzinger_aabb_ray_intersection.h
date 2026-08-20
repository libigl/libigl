// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_EYTZINGER_AABB_RAY_INTERSECTION_H
#define IGL_EYTZINGER_AABB_RAY_INTERSECTION_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <functional>
#include <vector>

namespace igl
{
  /// Traverse an Eytzinger AABB tree to find leaf primitives intersected by a
  /// ray. The per-primitive intersection test and the (exact) hit-ordering
  /// comparison are supplied as function handles so that this traversal carries
  /// no exact-valued state: the only thing propagated for the closest hit is the
  /// identity (leaf index) of the current best primitive. Box culling and
  /// front-to-back ordering use ordinary floating point and are made conservative
  /// (boxes are slightly inflated, distance pruning uses a safety margin) so that
  /// a genuine hit is never skipped.
  ///
  /// @tparam FirstOnly  if true, `hits` receives the single closest primitive (as
  ///   decided by `before`); otherwise `hits` receives every intersected
  ///   primitive sorted by `before`.
  /// @param[in] source  dim-vector origin of ray
  /// @param[in] dir     dim-vector direction of ray
  /// @param[in] hit     function handle `hit(i, t_approx)` returning true if the
  ///   ray intersects primitive i and, when true, setting `t_approx` (of the box
  ///   scalar type) to an approximate ray parameter used only for culling/ordering
  ///   (std::function so this routine can be instantiated with local lambdas)
  /// @param[in] before  function handle `before(i, j)` returning true if
  ///   primitive i's hit is strictly closer than primitive j's (exact)
  /// @param[in] B1  #B by dim list of minimum corners of the Eytzinger AABBs
  /// @param[in] B2  #B by dim list of maximum corners of the Eytzinger AABBs
  /// @param[in] leaf #B list of leaf indices, -1 indicates internal node, -2
  ///   indicates empty node
  /// @param[out] hits  list of intersected primitive indices (size ≤ 1 if
  ///   FirstOnly), sorted near-to-far by `before`
  ///
  /// \see eytzinger_aabb, exact_ray_mesh_intersect
  template <
    bool FirstOnly,
    typename Derivedsource,
    typename Deriveddir,
    typename DerivedB,
    typename Derivedleaf>
  IGL_INLINE void eytzinger_aabb_ray_intersection(
    const Eigen::MatrixBase<Derivedsource> & source,
    const Eigen::MatrixBase<Deriveddir> & dir,
    const std::function<bool(const int, typename DerivedB::Scalar &)> & hit,
    const std::function<bool(const int, const int)> & before,
    const Eigen::MatrixBase<DerivedB> & B1,
    const Eigen::MatrixBase<DerivedB> & B2,
    const Eigen::MatrixBase<Derivedleaf> & leaf,
    std::vector<int> & hits);
}

#ifndef IGL_STATIC_LIBRARY
#  include "eytzinger_aabb_ray_intersection.cpp"
#endif
#endif
