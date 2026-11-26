// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PREDICATES_POINT_MESH_CLOSEST_POINT_H
#define IGL_PREDICATES_POINT_MESH_CLOSEST_POINT_H
#include "../igl_inline.h"
#include <Eigen/Core>
namespace igl { namespace predicates {
  /// Find an exact closest feature and construct its approximate closest point.
  /// Uses one traversal of the supplied Eytzinger tree. Equal distances choose
  /// the lowest triangle id. Degenerate triangles are supported; feature ids
  /// follow exact_point_triangle_closest_feature (edge i joins i to (i+1)%3).
  /// Inside/outside classification and signed-distance wrappers additionally
  /// require nondegenerate faces, even though this unsigned search supports them.
  /// The feature is selected using exact predicates on the original inputs;
  /// numerical construction never reclassifies that feature.
  ///
  /// Distance, point coordinates and barycentrics are approximate doubles,
  /// not correctly rounded exact constructions. Barycentric coefficients are
  /// zero exactly at the inactive corners of the selected feature and strictly
  /// positive at its active corners. Thus one, two or three positive entries
  /// identify a vertex, open edge or open face, respectively. Positive weights
  /// that round to zero are replaced by the smallest positive double.
  /// Barycentrics need not sum to exactly one, and C need not lie exactly on
  /// the triangle or equal interpolation using the rounded barycentrics.
  /// Distance is evaluated from the selected feature's squared-distance
  /// fraction, independently of C, and may underflow/overflow. No expansions
  /// survive this call.
  ///
  /// Tree validity, input types and arithmetic assumptions are those of
  /// exact_point_mesh_closest_triangle. An empty mesh/tree throws
  /// std::invalid_argument because it has no closest point.
  /// @param[in] query  Query 3-vector
  /// @param[in] V,F  #V by 3 coordinates and #F by 3 triangle vertex indices
  /// @param[in] B1,B2,leaf  Precomputed conservative Eytzinger tree
  /// @param[out] I  Index of the selected closest triangle
  /// @param[out] C  Approximate closest point (3-vector)
  /// @param[out] bary  Barycentric 3-vector in F.row(I) corner order; its zero
  ///   pattern identifies the exact feature (edge i joins i to (i+1)%3)
  /// @return Approximate unsigned distance
  template<typename DerivedQ,typename DerivedV,typename DerivedF,
    typename DerivedB,typename DerivedL>
  IGL_INLINE double point_mesh_closest_point(
    const Eigen::MatrixBase<DerivedQ> &query,
    const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
    const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
    const Eigen::MatrixBase<DerivedL> &leaf,
    Eigen::Index &I,Eigen::RowVector3d &C,Eigen::RowVector3d &bary);
}}
#ifndef IGL_STATIC_LIBRARY
#  include "point_mesh_closest_point.cpp"
#endif
#endif
