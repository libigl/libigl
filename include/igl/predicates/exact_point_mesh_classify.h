// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PREDICATES_EXACT_POINT_MESH_CLASSIFY_H
#define IGL_PREDICATES_EXACT_POINT_MESH_CLASSIFY_H
#include "../igl_inline.h"
#include "exact_mesh_feature_exposure.h"
#include <Eigen/Core>
namespace igl { namespace predicates {
  enum class PointMeshLocation { INSIDE=-1, ON_SURFACE=0, OUTSIDE=1 };
  /// Classify a query using a triangle containing a globally closest point.
  /// The mesh must be embedded, closed, vertex-manifold, consistently oriented
  /// with normals pointing out of the solid, and have nondegenerate faces.
  /// Cavity shells must point into cavities. The triangle may be any exact
  /// minimizer, including ties. Supplying another triangle violates the contract.
  /// Uses orient3d on open faces, the adjacent fold on open edges, and the fold
  /// of an incident edge minimizing squared normalized projection at vertices.
  /// The latter choice uses metric predicates without normalizing any vector.
  /// Inputs and arithmetic follow exact_point_triangle_closest_feature's contract.
  /// @param[in] query  Query 3-vector
  /// @param[in] V,F  #V by 3 coordinates and #F by 3 vertex indices
  /// @param[in] closest_triangle  Exact minimizing triangle index
  /// @param[in] TT,TTi  Array output of triangle_triangle_adjacency(F,TT,TTi)
  /// @return INSIDE, OUTSIDE, or ON_SURFACE
  template<typename DerivedQ,typename DerivedV,typename DerivedF,
    typename DerivedTT,typename DerivedTTi>
  IGL_INLINE PointMeshLocation exact_point_mesh_classify(
    const Eigen::MatrixBase<DerivedQ> &query,
    const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
    Eigen::Index closest_triangle,const Eigen::MatrixBase<DerivedTT> &TT,
    const Eigen::MatrixBase<DerivedTTi> &TTi);
  /// Use a cache computed by exact_mesh_feature_exposure on the same mesh.
  /// @param[in] exposure  Current edge and vertex labels; replaces TT and TTi
  template<typename DerivedQ,typename DerivedV,typename DerivedF>
  IGL_INLINE PointMeshLocation exact_point_mesh_classify(
    const Eigen::MatrixBase<DerivedQ> &query,
    const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
    Eigen::Index closest_triangle,const MeshFeatureExposure &exposure);
}}
#ifndef IGL_STATIC_LIBRARY
#  include "exact_point_mesh_classify.cpp"
#endif
#endif
