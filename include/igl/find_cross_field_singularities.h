// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>, Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_FIND_CROSS_FIELD_SINGULARITIES_H
#define IGL_FIND_CROSS_FIELD_SINGULARITIES_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  /// Computes the Z/4 holonomy of a cross field around each vertex from its
  /// per-edge integer mismatch (see igl::cross_field_mismatch).
  ///
  /// `singularityIndex(v)` is the sum of the mismatches around the star of v
  /// reduced modulo 4 into {0,1,2,3}. This is the residue mod 4 of the true
  /// signed singularity index in quarter turns, which is exactly what seam
  /// computation (igl::cut_mesh_from_singularities) and integer-grid matching
  /// need, but it is not the index itself: -1 is reported as 3, -2 as 2, and a
  /// vertex around which the field makes a full turn (index 4 quarter turns,
  /// e.g. the pole of a longitude/latitude field on a sphere) is reported as 0
  /// and marked non-singular. Consequently the values do not sum to 4χ. For the
  /// signed integer index use igl::cross_field_singularity_index.
  ///
  /// @param[in] V                #V by 3 eigen Matrix of mesh vertex 3D positions
  /// @param[in] F                #F by 3 eigen Matrix of triangle indices
  /// @param[in] mismatch         #F by 3 eigen Matrix containing the integer mismatch of the cross field
  ///                    across all face edges
  /// @param[out] isSingularity    #V by 1 boolean eigen Vector indicating the presence of a singularity on a vertex
  ///                    (true iff singularityIndex is nonzero)
  /// @param[out] singularityIndex #V by 1 integer eigen Vector containing the singularity index mod 4
  ///                    in {0,1,2,3}. Boundary vertices are assigned 0.
  ///
  /// \see cross_field_singularity_index
  template <typename DerivedV, typename DerivedF, typename DerivedM, typename DerivedO>
  IGL_INLINE void find_cross_field_singularities(const Eigen::MatrixBase<DerivedV> &V,
                                                 const Eigen::MatrixBase<DerivedF> &F,
                                                 const Eigen::MatrixBase<DerivedM> &mismatch,
                                                 Eigen::PlainObjectBase<DerivedO> &isSingularity,
                                                 Eigen::PlainObjectBase<DerivedO> &singularityIndex);
  /// \overload
  ///
  /// \brief Wrapper that calculates the mismatch if it is not provided.
  ///
  /// @param[in] PD1              #F by 3 eigen Matrix of the first per face cross field vector
  /// @param[in] PD2              #F by 3 eigen Matrix of the second per face  cross field vector
  /// @param[in] isCombed        boolean indicating whether the cross field is already combed
  ///                    (see igl::comb_cross_field). If false, the field is combed internally.
  template <typename DerivedV, typename DerivedF, typename DerivedO>
  IGL_INLINE void find_cross_field_singularities(const Eigen::MatrixBase<DerivedV> &V,
                                                 const Eigen::MatrixBase<DerivedF> &F,
                                                 const Eigen::MatrixBase<DerivedV> &PD1,
                                                 const Eigen::MatrixBase<DerivedV> &PD2,
                                                 Eigen::PlainObjectBase<DerivedO> &isSingularity,
                                                 Eigen::PlainObjectBase<DerivedO> &singularityIndex,
                                                 bool isCombed = false);
}
#ifndef IGL_STATIC_LIBRARY
#include "find_cross_field_singularities.cpp"
#endif

#endif
