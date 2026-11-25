// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PREDICATES_EXACT_POINT_MESH_CLOSEST_TRIANGLE_H
#define IGL_PREDICATES_EXACT_POINT_MESH_CLOSEST_TRIANGLE_H
#include "../igl_inline.h"
#include <Eigen/Core>
namespace igl { namespace predicates {
  /// Find a triangle containing a globally closest point using an Eytzinger tree.
  /// Bounds must conservatively contain their triangles, and the tree must
  /// contain every face exactly once. Use igl::eytzinger_aabb to build it.
  /// Pruning uses outward-rounded distance bounds; leaf comparisons are exact.
  /// Degenerate triangles are supported. Equal distances choose the lowest id.
  /// Inside/outside classification and signed-distance wrappers additionally
  /// require nondegenerate faces, even though this unsigned search supports them.
  /// Inputs must be finite float or double coordinates. Link igl::predicates
  /// with its double backend. Correctness assumes round-to-nearest with gradual
  /// underflow, without reassociation or contraction. Exact refinement rescales
  /// coordinates by a common power of two and throws std::range_error if their
  /// exponent span cannot fit safely. A sufficient range for the entire point
  /// predicate family is: after rescaling, every nonzero coordinate's lowest
  /// set bit is at least 2^-100 and its magnitude is less than 2^91. Predicates
  /// of lower degree can refine wider spans. No expansion survives the call.
  /// @param[in] query  Query 3-vector
  /// @param[in] V,F  #V by 3 coordinates and #F by 3 triangle vertex indices
  /// @param[in] B1,B2,leaf  Eytzinger boxes and leaf ids (-1 internal, -2 empty)
  /// @return Closest face index, or -1 for an empty tree
  template<typename DerivedQ,typename DerivedV,typename DerivedF,
    typename DerivedB,typename DerivedL>
  IGL_INLINE Eigen::Index exact_point_mesh_closest_triangle(
    const Eigen::MatrixBase<DerivedQ> &query,
    const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
    const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
    const Eigen::MatrixBase<DerivedL> &leaf);
}}
#ifndef IGL_STATIC_LIBRARY
#  include "exact_point_mesh_closest_triangle.cpp"
#endif
#endif
