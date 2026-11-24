// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PREDICATES_EXACT_RAY_MESH_INTERSECT_H
#define IGL_PREDICATES_EXACT_RAY_MESH_INTERSECT_H
#include "../igl_inline.h"
#include "../Hit.h"
#include <Eigen/Core>
#include <vector>

namespace igl
{
namespace predicates
{
  /// Shoot a ray against a mesh (V,F) using exact leaf classification and hit
  /// ordering. The sign predicates refer to the mathematical ray defined by the
  /// floating point `source` and `dir`, without constructing an intersection
  /// location or ray parameter. Coordinates must be finite `float` or `double`.
  ///
  /// Candidate culling uses the floating-point eytzinger_aabb_ray_intersection
  /// routine and retains its numerical limitations: distant origins and very
  /// small direction components can cause missed candidates even when the
  /// leaf predicates would classify them correctly.
  ///
  /// Only clean transverse hits are reported. A triangle the ray is coplanar with
  /// (exact_ray_triangle_intersect returns Coplanar) may intersect in a segment.
  /// These predicates do not resolve coplanar intersections, so they are ignored
  /// (classified as no-hit).
  ///
  /// Each query builds an AABB tree. Runtime depends on the number of candidate
  /// triangles and hit comparisons, and on how often sign filters need exact
  /// expansion arithmetic.
  ///
  /// Link against igl::predicates with its double-precision backend. Correctness
  /// assumes round-to-nearest arithmetic without fast-math reassociation or
  /// intermediate overflow/underflow, as in Shewchuk's predicates. Float inputs
  /// are promoted exactly to double. No expansion is retained between calls.
  ///
  /// @param[in] source  3-vector origin of ray
  /// @param[in] dir     3-vector direction of ray
  /// @param[in] V  #V by 3 list of mesh vertex positions
  /// @param[in] F  #F by 3 list of mesh face indices into V
  /// @param[out] fid  index into F of the first hit face (set only if returning true)
  /// @return true if there was a hit
  ///
  /// \see exact_ray_triangle_intersect, igl::ray_mesh_intersect
  template <
    typename Derivedsource,
    typename Deriveddir,
    typename DerivedV,
    typename DerivedF>
  IGL_INLINE bool exact_ray_mesh_intersect(
    const Eigen::MatrixBase<Derivedsource> & source,
    const Eigen::MatrixBase<Deriveddir> & dir,
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    int & fid);
  /// Collect every triangle transversally intersected by the ray, sorted exactly
  /// near-to-far.
  ///
  /// @param[out] fids  indices into F of hit faces, sorted by exact hit order
  /// @return true if there were any hits
  template <
    typename Derivedsource,
    typename Deriveddir,
    typename DerivedV,
    typename DerivedF>
  IGL_INLINE bool exact_ray_mesh_intersect_all(
    const Eigen::MatrixBase<Derivedsource> & source,
    const Eigen::MatrixBase<Deriveddir> & dir,
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    std::vector<int> & fids);
  /// \overload
  /// Convenience overload that additionally constructs an approximate (floating
  /// point) hit record for the first hit as a post-process. The returned
  /// barycentric coordinates and parameter t are computed inexactly and can be
  /// unreliable on ill-conditioned triangles, including zero values when the
  /// rounded determinant vanishes.
  ///
  /// @param[out] hit  first hit (id, u, v, t), set only if returning true
  template <
    typename Derivedsource,
    typename Deriveddir,
    typename DerivedV,
    typename DerivedF>
  IGL_INLINE bool exact_ray_mesh_intersect(
    const Eigen::MatrixBase<Derivedsource> & source,
    const Eigen::MatrixBase<Deriveddir> & dir,
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    igl::Hit<typename DerivedV::Scalar> & hit);
  /// \overload
  /// Convenience overload constructing approximate hit records for all hits,
  /// sorted by exact hit order.
  ///
  /// @param[out] hits  sorted list of hits with approximate (u,v,t)
  template <
    typename Derivedsource,
    typename Deriveddir,
    typename DerivedV,
    typename DerivedF>
  IGL_INLINE bool exact_ray_mesh_intersect_all(
    const Eigen::MatrixBase<Derivedsource> & source,
    const Eigen::MatrixBase<Deriveddir> & dir,
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    std::vector<igl::Hit<typename DerivedV::Scalar>> & hits);
}
}

#ifndef IGL_STATIC_LIBRARY
#  include "exact_ray_mesh_intersect.cpp"
#endif
#endif
