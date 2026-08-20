// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_EXACT_RAY_MESH_INTERSECT_H
#define IGL_EXACT_RAY_MESH_INTERSECT_H
#include "igl_inline.h"
#include "Hit.h"
#include <Eigen/Core>
#include <vector>

namespace igl
{
  /// Shoot a ray against a mesh (V,F) and robustly determine the first (closest)
  /// intersected triangle using only exact sign predicates. The classification —
  /// which triangle is hit and in what order — is exact with respect to the
  /// mathematical ray defined by the floating point `source` and `dir`; there are
  /// no tolerances or fudge factors, and no intersection location or ray parameter
  /// is ever constructed for that decision (constructing it would be inexact).
  /// Coordinates must be `float` or `double`.
  ///
  /// Only clean transverse hits are reported. A triangle the ray is coplanar with
  /// (exact_ray_triangle_intersect returns Coplanar) has a segment-valued
  /// intersection that cannot be ordered against point hits with these predicates,
  /// so it is deliberately ignored (classified as no-hit).
  ///
  /// This exactness is expensive: the exact per-triangle predicates dominate, so
  /// as a rough guide it runs one to two orders of magnitude slower than
  /// igl::ray_mesh_intersect (~40× on a 1k-face mesh in local benchmarks; the
  /// AABB tree recovers some of that for first-hit queries on larger meshes).
  /// Prefer igl::ray_mesh_intersect / igl::AABB unless exact classification is
  /// required.
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
  /// point) hit record for the first hit as a post-process. The classification
  /// (which face, and that it is the closest) is exact; the returned barycentric
  /// coordinates and parameter t are computed inexactly for convenience only.
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

#ifndef IGL_STATIC_LIBRARY
#  include "exact_ray_mesh_intersect.cpp"
#endif
#endif
