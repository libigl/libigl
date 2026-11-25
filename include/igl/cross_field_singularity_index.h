// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Derek Wingard
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_CROSS_FIELD_SINGULARITY_INDEX_H
#define IGL_CROSS_FIELD_SINGULARITY_INDEX_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  /// Signed integer singularity index of a per-face N-RoSy (rotationally
  /// symmetric) direction field at every vertex, in units of 1/N turns. For a
  /// cross field (N = 4) the unit is a quarter turn, so a vertex of index +1/4
  /// (the typical valence-3 quad-mesh singularity) reports 1, a vertex of index
  /// -1/4 reports -1, and a vertex around which the field makes a full turn
  /// (e.g. the pole of a longitude/latitude field on a sphere) reports 4.
  ///
  /// The index at vertex v is
  ///
  ///     index(v) = N * (K_v + Σ_e r_e) / (2π)
  ///
  /// where K_v is the angle defect at v and r_e ∈ [-π/N, π/N) is the residual
  /// rotation of the field across each edge e of the star of v after
  /// transporting it across the edge (unfolding) and removing the nearest
  /// multiple of 2π/N. The quantity on the right is an exact integer for every
  /// field; it is rounded only to absorb floating point error. Because r is
  /// antisymmetric and Σ_v K_v = 2πχ (discrete Gauss–Bonnet), the indices of a
  /// closed mesh always satisfy
  ///
  ///     Σ_v index(v) = N χ(M)
  ///
  /// (Poincaré–Hopf). This differs from igl::find_cross_field_singularities,
  /// whose output is only the residue index(v) mod 4 stored in {0,1,2,3}: that
  /// representation cannot distinguish -1 from 3, 2 from -2, or 4 from 0.
  ///
  /// Only one representative direction per face is needed; the field does not
  /// have to be combed, since the residual r_e does not depend on which of the
  /// N symmetric directions is stored.
  ///
  /// The index is a faithful description of the field only when the mesh
  /// resolves it, i.e. when the field turns by less than π/N relative to
  /// parallel transport across every edge of the star. A field that turns
  /// faster than that is aliased (as it is for any discrete method), but the
  /// sum rule above still holds.
  ///
  /// @param[in] V    #V by 3 list of vertex positions
  /// @param[in] F    #F by 3 list of triangle indices, consistently
  ///                 oriented, edge-manifold
  /// @param[in] PD1  #F by 3 list of one representative direction of the
  ///                 field per face (need not be unit length or exactly
  ///                 tangent; it is projected onto the face)
  /// @param[in] N    rotational symmetry degree of the field: 4 for a cross
  ///                 field, 2 for a line field, 1 for a vector field
  /// @param[out] index  #V list of signed integer indices in units of 1/N
  ///                 turns. Boundary vertices are assigned 0.
  ///
  /// \see find_cross_field_singularities, cross_field_mismatch
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedP,
    typename DerivedI>
  IGL_INLINE void cross_field_singularity_index(
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    const Eigen::MatrixBase<DerivedP> & PD1,
    const int N,
    Eigen::PlainObjectBase<DerivedI> & index);
  /// \overload
  /// \brief Cross field (N = 4): index in quarter turns.
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedP,
    typename DerivedI>
  IGL_INLINE void cross_field_singularity_index(
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    const Eigen::MatrixBase<DerivedP> & PD1,
    Eigen::PlainObjectBase<DerivedI> & index);
}

#ifndef IGL_STATIC_LIBRARY
#  include "cross_field_singularity_index.cpp"
#endif

#endif
