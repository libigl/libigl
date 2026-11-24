// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PREDICATES_APPROXIMATE_SIGNED_DISTANCE_EXACT_SIGN_H
#define IGL_PREDICATES_APPROXIMATE_SIGNED_DISTANCE_EXACT_SIGN_H
#include "../igl_inline.h"
#include "point_mesh_closest_point.h"
#include "exact_point_mesh_classify.h"
#include <vector>
namespace igl { namespace predicates {
  /// Compute an approximate signed distance with exact inside/outside location.
  /// Performs one closest-point traversal and signs its selected feature using
  /// the supplied exposure cache. Inside distances are negative; outside
  /// distances are positive. ON_SURFACE returns exactly +0.0. The location
  /// output remains authoritative if a nonzero magnitude rounds to zero.
  ///
  /// Requires a nonempty, embedded, closed, vertex-manifold triangle mesh with
  /// nondegenerate faces oriented out of the solid. Build exposure with
  /// exact_mesh_feature_exposure using array triangle_triangle_adjacency.
  /// Rebuild both tree and exposure after changing geometry/topology.
  /// Global mesh validity and cache freshness are the caller's responsibility.
  /// Exactness has the input/range/rounding assumptions documented by
  /// exact_point_mesh_closest_triangle and exact_point_mesh_classify.
  /// The magnitude and numerical closest-point data are approximate, as
  /// documented by point_mesh_closest_point. No expansions are retained.
  /// @param[in] query  Query 3-vector
  /// @param[in] V,F  #V by 3 coordinates and #F by 3 triangle vertex indices
  /// @param[in] B1,B2,leaf  Precomputed conservative Eytzinger tree
  /// @param[in] exposure  Precomputed feature exposure cache for V,F
  /// @param[out] closest  Exact triangle/feature ids and approximate point,
  ///   barycentrics, unsigned distance and squared distance
  /// @param[out] location  Exact INSIDE, ON_SURFACE or OUTSIDE classification
  /// @return Approximate signed distance
  template<typename DerivedQ,typename DerivedV,typename DerivedF,
    typename DerivedB,typename DerivedL>
  IGL_INLINE double approximate_signed_distance_exact_sign(
    const Eigen::MatrixBase<DerivedQ> &query,
    const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
    const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
    const Eigen::MatrixBase<DerivedL> &leaf,const MeshFeatureExposure &exposure,
    PointMeshClosestPoint &closest,PointMeshLocation &location);

  /// Batched overload reusing the tree and exposure cache for every query.
  /// The single-query contract applies to each row. Outputs are resized;
  /// an empty query matrix gives empty outputs. An empty mesh/tree throws
  /// std::invalid_argument. Output arrays must not alias input storage.
  /// @param[in] P  #P by 3 query positions
  /// @param[in] V,F,B1,B2,leaf,exposure  As in the single-query overload
  /// @param[out] S  #P approximate signed distances
  /// @param[out] closest  #P closest-point results, corresponding to rows of P
  /// @param[out] locations  #P exact location values (-1 inside, 0 surface,
  ///   +1 outside), independent of rounding in S
  template<typename DerivedP,typename DerivedV,typename DerivedF,
    typename DerivedB,typename DerivedL>
  IGL_INLINE void approximate_signed_distance_exact_sign(
    const Eigen::MatrixBase<DerivedP> &P,
    const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
    const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
    const Eigen::MatrixBase<DerivedL> &leaf,const MeshFeatureExposure &exposure,
    Eigen::VectorXd &S,std::vector<PointMeshClosestPoint> &closest,
    Eigen::VectorXi &locations);
}}
#ifndef IGL_STATIC_LIBRARY
#  include "approximate_signed_distance_exact_sign.cpp"
#endif
#endif
