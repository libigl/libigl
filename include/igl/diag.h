// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_DIAG_H
#define IGL_DIAG_H
#include "igl_inline.h"
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>

namespace igl
{
  /// Extracts the main diagonal of a matrix as a vector. Like matlab's `diag`.
  ///
  ///
  /// @tparam T  should be a eigen sparse matrix primitive type like int or double
  /// @param[in] X  an m by n sparse matrix
  /// @param[out] V  a min(m,n) sparse vector
  ///
  /// http://forum.kde.org/viewtopic.php?f=74&t=117476&p=292388#p292388
  ///
  /// \deprecated Use 
  ///   `VectorXd V = X.diagonal()` and 
  ///   `SparseVector<double> V = X.diagonal().sparseView()`
  ///   `SparseMatrix<double> X = V.asDiagonal().sparseView()`
  ///
  ///
  template <typename T>
  IGL_INLINE void diag(
    const Eigen::SparseMatrix<T>& X, 
    Eigen::SparseVector<T>& V);
  /// \overload
  template <typename T,typename DerivedV>
  IGL_INLINE void diag(
    const Eigen::SparseMatrix<T>& X, 
    Eigen::MatrixBase<DerivedV>& V);
  /// Builds a sparse matrix with a given vector along the main diagonal.
  ///
  /// @tparam T  should be a eigen sparse matrix primitive type like int or double
  /// @param[in] V  a m sparse vector
  /// @param[out] X  a m by m sparse matrix
  template <typename T>
  IGL_INLINE void diag(
    const Eigen::SparseVector<T>& V,
    Eigen::SparseMatrix<T>& X);
  /// \overload
  template <typename T, typename DerivedV>
  IGL_INLINE void diag(
    const Eigen::MatrixBase<DerivedV>& V,
    Eigen::SparseMatrix<T>& X);
}

#ifndef IGL_STATIC_LIBRARY
#  include "diag.cpp"
#endif

#endif
