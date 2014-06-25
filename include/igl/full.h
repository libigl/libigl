// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_FULL_H
#define IGL_FULL_H
#include "igl_inline.h"
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>
namespace igl
{
  // This is totally unnecessary. You can just call MatrixXd B = MatrixXd(A);
  //
  // Convert a sparsematrix into a full one
  //
  // Templates:
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Input:
  //   A  m by n sparse matrix
  // Output:
  //   B  m by n dense/full matrix
  template <typename T,typename DerivedB>
  IGL_INLINE void full(
    const Eigen::SparseMatrix<T> & A,
    Eigen::PlainObjectBase<DerivedB> & B);
  // If already full then this will just be a copy by assignment
  template <typename DerivedA,typename DerivedB>
  IGL_INLINE void full(
    const Eigen::PlainObjectBase<DerivedA>& A,
    Eigen::PlainObjectBase<DerivedB>& B);
}

#ifndef IGL_STATIC_LIBRARY
#  include "full.cpp"
#endif

#endif
