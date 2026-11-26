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
namespace igl { namespace predicates {
  /// Compute an approximate signed distance with exact inside/outside location.
  /// Performs one closest-point traversal and signs its selected feature using
  /// the supplied exposure cache. Inside distances are negative; outside
  /// distances are positive; exactly on-surface queries return +0.0. A nonzero
  /// magnitude that rounds to zero is replaced by the smallest positive double
  /// before applying the sign. Thus the sign and zeros of S retain the exact
  /// inside/on-surface/outside classification, even after underflow.
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
  /// @param[out] I  Index of the selected closest triangle; exact ties choose
  ///   the lowest triangle id
  /// @param[out] C  Approximate closest point (3-vector)
  /// @param[out] bary  Barycentric 3-vector in F.row(I) corner order. Its zero
  ///   pattern identifies the exact feature as in point_mesh_closest_point.
  /// @return Approximate signed distance
  template<typename DerivedQ,typename DerivedV,typename DerivedF,
    typename DerivedB,typename DerivedL>
  IGL_INLINE double approximate_signed_distance_exact_sign(
    const Eigen::MatrixBase<DerivedQ> &query,
    const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
    const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
    const Eigen::MatrixBase<DerivedL> &leaf,const MeshFeatureExposure &exposure,
    Eigen::Index &I,Eigen::RowVector3d &C,Eigen::RowVector3d &bary);

  /// Uncached overload. Reuse array triangle_triangle_adjacency outputs TT,TTi
  /// instead of an exposure cache; evaluate local edge/vertex signs per query.
  /// Requires nondegenerate faces, as well as the cached overload's other mesh,
  /// numerical and output contracts. Degenerate incident faces encountered at
  /// a zero fold throw std::invalid_argument; this is not a full mesh validation.
  /// @param[in] TT,TTi  #F by 3 neighboring faces and corresponding edge ids
  template<typename DerivedQ,typename DerivedV,typename DerivedF,
    typename DerivedB,typename DerivedL,typename DerivedTT,typename DerivedTTi>
  IGL_INLINE double approximate_signed_distance_exact_sign(
    const Eigen::MatrixBase<DerivedQ> &query,
    const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
    const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
    const Eigen::MatrixBase<DerivedL> &leaf,const Eigen::MatrixBase<DerivedTT> &TT,
    const Eigen::MatrixBase<DerivedTTi> &TTi,
    Eigen::Index &I,Eigen::RowVector3d &C,Eigen::RowVector3d &bary);

  /// Batched overload reusing the tree and exposure cache for every query.
  /// The single-query contract applies to each row. Outputs are resized;
  /// an empty query matrix gives empty outputs. An empty mesh/tree throws
  /// std::invalid_argument. Output arrays must not alias each other or input storage.
  /// Batches of at least 10,000 queries use igl::parallel_for with the configured
  /// igl::default_num_threads(). Worker exceptions are rethrown on the caller
  /// after joining all workers; outputs may be partially filled on failure.
  /// @param[in] P  #P by 3 query positions
  /// @param[in] V,F,B1,B2,leaf,exposure  As in the single-query overload
  /// @param[out] S  #P approximate signed distances
  /// @param[out] I  #P selected closest triangle indices
  /// @param[out] C  #P by 3 approximate closest points
  /// @param[out] bary  #P by 3 barycentrics; each row's zero pattern identifies
  ///   its exact closest feature
  template<typename DerivedP,typename DerivedV,typename DerivedF,
    typename DerivedB,typename DerivedL>
  IGL_INLINE void approximate_signed_distance_exact_sign(
    const Eigen::MatrixBase<DerivedP> &P,
    const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
    const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
    const Eigen::MatrixBase<DerivedL> &leaf,const MeshFeatureExposure &exposure,
    Eigen::VectorXd &S,Eigen::VectorXi &I,Eigen::MatrixXd &C,Eigen::MatrixXd &bary);

  /// Uncached batched overload with the same parallel policy and outputs.
  template<typename DerivedP,typename DerivedV,typename DerivedF,
    typename DerivedB,typename DerivedL,typename DerivedTT,typename DerivedTTi>
  IGL_INLINE void approximate_signed_distance_exact_sign(
    const Eigen::MatrixBase<DerivedP> &P,
    const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
    const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
    const Eigen::MatrixBase<DerivedL> &leaf,const Eigen::MatrixBase<DerivedTT> &TT,
    const Eigen::MatrixBase<DerivedTTi> &TTi,Eigen::VectorXd &S,
    Eigen::VectorXi &I,Eigen::MatrixXd &C,Eigen::MatrixXd &bary);
}}
#ifndef IGL_STATIC_LIBRARY
#  include "approximate_signed_distance_exact_sign.cpp"
#endif
#endif
