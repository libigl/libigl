// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2020 Oded Stein <oded.stein@columbia.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_CR_VECTOR_LAPLACIAN_H
#define IGL_CR_VECTOR_LAPLACIAN_H

#include "igl_inline.h"

#include <Eigen/Core>
#include <Eigen/Sparse>


namespace igl
{
  /// Computes the CR vector Laplacian matrix.
  /// See Oded Stein, Max Wardetzky, Alec Jacobson, Eitan Grinspun, 2020.
  ///  "A Simple Discretization of the Vector Dirichlet Energy"
  ///
  ///  @param[in] V #V by 3 list of mesh vertex positions
  ///  @param[in] F #F by 3 list of mesh face indices into rows of V
  ///  @param[in] E #F by 3 a mapping from each halfedge to each edge
  ///  @param[in] oE #F by 3 the orientation (e.g., -1 or 1) of each halfedge
  ///    compared to the orientation of the actual edge, as computed with
  ///    orient_halfedges. will be computed if not provided.
  ///  @param[out] L 2*|HE| by 2*|HE| computed Laplacian matrix
  template <typename DerivedV, typename DerivedF, typename DerivedE,
  typename DerivedOE, typename ScalarL>
  IGL_INLINE void
  cr_vector_laplacian(
    const Eigen::MatrixBase<DerivedV>& V,
    const Eigen::MatrixBase<DerivedF>& F,
    const Eigen::MatrixBase<DerivedE>& E,
    const Eigen::MatrixBase<DerivedOE>& oE,
    Eigen::SparseMatrix<ScalarL>& L);
  /// \overload
  ///
  /// \brief `E` and `oE` are computed and output.
  template <typename DerivedV, typename DerivedF, typename DerivedE,
  typename DerivedOE, typename ScalarL>
  IGL_INLINE void
  cr_vector_laplacian(
    const Eigen::MatrixBase<DerivedV>& V,
    const Eigen::MatrixBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedE>& E,
    Eigen::PlainObjectBase<DerivedOE>& oE,
    Eigen::SparseMatrix<ScalarL>& L);
  /// \overload
  /// \brief intrinsic version.
  ///
  ///  @param[in] l_sq #F by 3 list of squared edge lengths of each halfedge
  ///  @param[in] dA #F list of double areas
  ///
  ///  \fileinfo
  template <typename DerivedF, typename DerivedL_sq, typename DeriveddA,
  typename DerivedE, typename DerivedOE, typename ScalarL>
  IGL_INLINE void
  cr_vector_laplacian_intrinsic(
    const Eigen::MatrixBase<DerivedF>& F,
    const Eigen::MatrixBase<DerivedL_sq>& l_sq,
    const Eigen::MatrixBase<DeriveddA>& dA,
    const Eigen::MatrixBase<DerivedE>& E,
    const Eigen::MatrixBase<DerivedOE>& oE,
    Eigen::SparseMatrix<ScalarL>& L);
  /// \overload
  /// \fileinfo
  template <typename DerivedF, typename DerivedL_sq, typename DerivedE,
  typename DerivedOE, typename ScalarL>
  IGL_INLINE void
  cr_vector_laplacian_intrinsic(
    const Eigen::MatrixBase<DerivedF>& F,
    const Eigen::MatrixBase<DerivedL_sq>& l_sq,
    const Eigen::MatrixBase<DerivedE>& E,
    const Eigen::MatrixBase<DerivedOE>& oE,
    Eigen::SparseMatrix<ScalarL>& L);
}


#ifndef IGL_STATIC_LIBRARY
#  include "cr_vector_laplacian.cpp"
#endif

#endif
