// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PREDICATES_POINT_MESH_CLOSEST_POINT_H
#define IGL_PREDICATES_POINT_MESH_CLOSEST_POINT_H
#include "../igl_inline.h"
#include "exact_point_triangle_closest_feature.h"
#include <Eigen/Core>
namespace igl { namespace predicates {
  /// Exact feature selection with approximate floating-point constructions.
  struct PointMeshClosestPoint
  {
    Eigen::Index triangle;
    ClosestPointFeature feature;
    int feature_id;
    double distance;
    double squared_distance;
    Eigen::RowVector3d closest_point;
    Eigen::RowVector3d barycentric_coordinates;
  };

  /// Find an exact closest feature and construct its approximate closest point.
  /// Uses one traversal of the supplied Eytzinger tree. Equal distances choose
  /// the lowest triangle id. Degenerate triangles are supported; feature ids
  /// follow exact_point_triangle_closest_feature (edge i joins i to (i+1)%3).
  /// The feature is selected using exact predicates on the original inputs;
  /// numerical construction never reclassifies that feature.
  ///
  /// Distances, point coordinates and barycentrics are approximate doubles,
  /// not correctly rounded exact constructions. Rounded barycentrics need not
  /// sum to exactly one or retain an open feature's strict inequalities, and
  /// closest_point need not lie exactly on the triangle. Distance is evaluated
  /// from the selected feature's squared-distance fraction, independently of
  /// the rounded closest point. squared_distance may underflow/overflow even
  /// when distance is representable. No expansions survive this call.
  ///
  /// Tree validity, input types and arithmetic assumptions are those of
  /// exact_point_mesh_closest_triangle. An empty mesh/tree throws
  /// std::invalid_argument because it has no closest point.
  /// @param[in] query  Query 3-vector
  /// @param[in] V,F  #V by 3 coordinates and #F by 3 triangle vertex indices
  /// @param[in] B1,B2,leaf  Precomputed conservative Eytzinger tree
  /// @return Exact triangle/feature ids and approximate numerical results
  template<typename DerivedQ,typename DerivedV,typename DerivedF,
    typename DerivedB,typename DerivedL>
  IGL_INLINE PointMeshClosestPoint point_mesh_closest_point(
    const Eigen::MatrixBase<DerivedQ> &query,
    const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
    const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
    const Eigen::MatrixBase<DerivedL> &leaf);
}}
#ifndef IGL_STATIC_LIBRARY
#  include "point_mesh_closest_point.cpp"
#endif
#endif
