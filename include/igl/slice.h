// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SLICE_H
#define IGL_SLICE_H
#include "igl_inline.h"

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>

namespace igl
{
  // THIS MAY HAVE BEEN SUPERSEDED BY EIGEN'S select FUNCTION
  // 
  // Act like the matlab X(row_indices,col_indices) operator
  // 
  // Inputs:
  //   X  m by n matrix
  //   R  list of row indices
  //   C  list of column indices
  // Output:
  //   Y  #R by #C matrix
  template <typename T>
  IGL_INLINE void slice(
    const Eigen::SparseMatrix<T>& X,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & C,
    Eigen::SparseMatrix<T>& Y);
  // Wrapper to only slice in one direction
  //
  // Inputs:
  //   dim  dimension to slice in 1 or 2, dim=1 --> X(R,:), dim=2 --> X(:,R)
  //
  // Note: For now this is just a cheap wrapper.
  template <typename Mat>
  IGL_INLINE void slice(
    const Mat& X,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
    const int dim,
    Mat& Y);
  template <typename DerivedX>
  IGL_INLINE void slice(
    const Eigen::PlainObjectBase<DerivedX> & X,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & C,
    Eigen::PlainObjectBase<DerivedX> & Y);
  template <typename DerivedX>
  IGL_INLINE void slice(
    const Eigen::PlainObjectBase<DerivedX> & X,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
    Eigen::PlainObjectBase<DerivedX> & Y);
  // VectorXi Y = slice(X,R);
  template <typename DerivedX>
  IGL_INLINE Eigen::PlainObjectBase<DerivedX> slice(
    const Eigen::PlainObjectBase<DerivedX> & X,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & R);
  template <typename DerivedX>
  IGL_INLINE Eigen::PlainObjectBase<DerivedX> slice(
    const Eigen::PlainObjectBase<DerivedX>& X,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
    const int dim);
}

#ifndef IGL_STATIC_LIBRARY
#  include "slice.cpp"
#endif

#endif
