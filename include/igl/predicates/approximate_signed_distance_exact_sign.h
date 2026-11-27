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
  /// Compute a signed distance with exact closest-feature selection and sign.
  /// Closest-point selection is exact: I identifies a triangle containing a
  /// point that minimizes distance to the query over the entire mesh. Within
  /// that triangle, the zero pattern of bary identifies the vertex, open edge
  /// or open face containing this exact closest point. These cases have one,
  /// two or three strictly positive barycentric entries, respectively; all
  /// other entries are exactly zero. Active weights that round to zero are
  /// replaced by the smallest positive value of the output scalar type to
  /// preserve this classification.
  ///
  /// The numerical barycentric values and distance magnitude are approximate.
  /// Reconstructing coordinates from I and bary therefore gives an approximate
  /// closest point, while the selected triangle and feature are exact. See
  /// point_mesh_closest_point for the numerical construction contract.
  ///
  /// Performs one closest-point traversal and signs the selected feature using
  /// the supplied exposure cache. Inside distances are negative; outside
  /// distances are positive; exactly on-surface queries return +0.0. A nonzero
  /// magnitude that rounds to zero is replaced by the smallest positive value
  /// of the distance output type before applying the sign. Thus the sign and
  /// zeros of S retain the exact inside/on-surface/outside classification,
  /// even after underflow.
  ///
  /// Requires a nonempty, embedded, closed, vertex-manifold triangle mesh with
  /// nondegenerate faces oriented out of the solid. Build exposure with
  /// exact_mesh_feature_exposure using array triangle_triangle_adjacency.
  /// Rebuild both tree and exposure after changing geometry/topology.
  /// Global mesh validity and cache freshness are the caller's responsibility.
  /// Exactness has the input/range/rounding assumptions documented by
  /// exact_point_mesh_closest_triangle and exact_point_mesh_classify.
  /// No expansions are retained.
  /// @param[in] query  Query 3-vector
  /// @param[in] V,F  #V by 3 coordinates and #F by 3 triangle vertex indices
  /// @param[in] B1,B2,leaf  Precomputed conservative Eytzinger tree
  /// @param[in] exposure  Precomputed feature exposure cache for V,F
  /// @param[out] I  Index of a triangle containing an exact closest point on
  ///   the mesh; exact ties choose the lowest triangle id
  /// @param[out] bary  Barycentric row or column 3-vector (resized) in F.row(I)
  ///   corner order. Its zero pattern identifies the feature containing the
  ///   exact closest point; its positive numerical values are approximate.
  /// @return Approximate signed distance in double precision
  template<typename DerivedQ,typename DerivedV,typename DerivedF,
    typename DerivedB,typename DerivedL,typename DerivedBary>
  IGL_INLINE double approximate_signed_distance_exact_sign(
    const Eigen::MatrixBase<DerivedQ> &query,
    const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
    const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
    const Eigen::MatrixBase<DerivedL> &leaf,const MeshFeatureExposure &exposure,
    Eigen::Index &I,Eigen::PlainObjectBase<DerivedBary> &bary);

  /// Uncached overload. Reuse array triangle_triangle_adjacency outputs TT,TTi
  /// instead of an exposure cache; evaluate local edge/vertex signs per query.
  /// Requires nondegenerate faces, as well as the cached overload's other mesh,
  /// numerical and output contracts. Degenerate incident faces encountered at
  /// a zero fold throw std::invalid_argument; this is not a full mesh validation.
  /// @param[in] TT,TTi  #F by 3 neighboring faces and corresponding edge ids
  template<typename DerivedQ,typename DerivedV,typename DerivedF,
    typename DerivedB,typename DerivedL,typename DerivedTT,typename DerivedTTi,typename DerivedBary>
  IGL_INLINE double approximate_signed_distance_exact_sign(
    const Eigen::MatrixBase<DerivedQ> &query,
    const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
    const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
    const Eigen::MatrixBase<DerivedL> &leaf,const Eigen::MatrixBase<DerivedTT> &TT,
    const Eigen::MatrixBase<DerivedTTi> &TTi,
    Eigen::Index &I,Eigen::PlainObjectBase<DerivedBary> &bary);

  /// Batched overload reusing the tree and exposure cache for every query.
  /// The single-query contract applies to each row. Outputs are resized;
  /// an empty query matrix gives empty outputs. An empty mesh/tree throws
  /// std::invalid_argument. Output arrays must not alias each other or input storage.
  /// Batches of at least 10,000 queries use igl::parallel_for with the configured
  /// igl::default_num_threads(). Worker exceptions are rethrown on the caller
  /// after joining all workers; outputs may be partially filled on failure.
  /// @param[in] P  #P by 3 query positions
  /// @param[in] V,F,B1,B2,leaf,exposure  As in the single-query overload
  /// @param[out] S  #P by 1 approximate signed distances; floating-point scalar
  /// @param[out] I  #P by 1 indices of triangles containing exact closest
  ///   points on the mesh; integral scalar capable of representing all indices
  /// @param[out] bary  #P by 3 approximate barycentrics; row i's zero pattern
  ///   identifies the feature of triangle I(i) containing an exact closest
  ///   point to P.row(i). Uses a floating-point scalar.
  template<typename DerivedP,typename DerivedV,typename DerivedF,
    typename DerivedB,typename DerivedL,
    typename DerivedS,typename DerivedI,typename DerivedBary>
  IGL_INLINE void approximate_signed_distance_exact_sign(
    const Eigen::MatrixBase<DerivedP> &P,
    const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
    const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
    const Eigen::MatrixBase<DerivedL> &leaf,const MeshFeatureExposure &exposure,
    Eigen::PlainObjectBase<DerivedS> &S,Eigen::PlainObjectBase<DerivedI> &I,
    Eigen::PlainObjectBase<DerivedBary> &bary);

  /// Uncached batched overload with the same parallel policy and outputs.
  template<typename DerivedP,typename DerivedV,typename DerivedF,
    typename DerivedB,typename DerivedL,typename DerivedTT,typename DerivedTTi,
    typename DerivedS,typename DerivedI,typename DerivedBary>
  IGL_INLINE void approximate_signed_distance_exact_sign(
    const Eigen::MatrixBase<DerivedP> &P,
    const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
    const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
    const Eigen::MatrixBase<DerivedL> &leaf,const Eigen::MatrixBase<DerivedTT> &TT,
    const Eigen::MatrixBase<DerivedTTi> &TTi,Eigen::PlainObjectBase<DerivedS> &S,
    Eigen::PlainObjectBase<DerivedI> &I,Eigen::PlainObjectBase<DerivedBary> &bary);
}}
#ifndef IGL_STATIC_LIBRARY
#  include "approximate_signed_distance_exact_sign.cpp"
#endif
#endif
