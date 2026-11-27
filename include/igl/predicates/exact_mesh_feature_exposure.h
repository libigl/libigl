// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PREDICATES_EXACT_MESH_FEATURE_EXPOSURE_H
#define IGL_PREDICATES_EXACT_MESH_FEATURE_EXPOSURE_H
#include "../igl_inline.h"
#include <Eigen/Core>
namespace igl { namespace predicates {
  enum class FeatureExposure
  {
    EXPOSED_TO_INSIDE=-1, FLAT=0, EXPOSED_TO_OUTSIDE=1, UNREACHABLE=2
  };
  /// Query-independent exposure of the closest-point features of a mesh.
  /// Edge (f,i) joins F(f,i) to F(f,(i+1)%3), matching array adjacency.
  /// Entries hold FeatureExposure values. Unreferenced vertices are UNREACHABLE.
  /// Invalidate this cache whenever coordinates, topology, or orientation change.
  struct MeshFeatureExposure
  {
    Eigen::MatrixXi edges;
    Eigen::VectorXi vertices;
  };
  /// Cache which side can reach each open edge and vertex as a closest point.
  /// The mesh must be embedded, closed, consistently outward-oriented, and
  /// vertex-manifold, with nondegenerate triangular faces. Connectivity is the
  /// array output of triangle_triangle_adjacency(F,TT,TTi).
  /// FLAT vertices have planar one-rings. UNREACHABLE vertices are strictly
  /// inside the 3D convex hull of their neighbors and themselves; vertices on
  /// its boundary are reachable. All tests use exact orientation signs.
  /// The implementation enumerates supporting planes in O(sum_v degree(v)^3).
  /// Arithmetic has the same range and rounding contract as Shewchuk orient3d.
  /// @param[in] V,F  #V by 3 coordinates and #F by 3 triangle vertex indices
  /// @param[in] TT,TTi  #F by 3 adjacent face and matching directed-edge ids
  /// @param[out] exposure  Edge and vertex exposure arrays
  template<typename DerivedV,typename DerivedF,typename DerivedTT,typename DerivedTTi>
  IGL_INLINE void exact_mesh_feature_exposure(
    const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
    const Eigen::MatrixBase<DerivedTT> &TT,const Eigen::MatrixBase<DerivedTTi> &TTi,
    MeshFeatureExposure &exposure);
}}
#ifndef IGL_STATIC_LIBRARY
#  include "exact_mesh_feature_exposure.cpp"
#endif
#endif
