// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_LU_LAGRANGE_H
#define IGL_LU_LAGRANGE_H
#include "igl_inline.h"

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>
namespace igl
{
  // KNOWN BUGS: This does not seem to be correct for non-empty C
  //
  // LU_LAGRANGE Compute a LU decomposition for a special type of
  // matrix Q that is symmetric but not positive-definite:
  // Q = [A'*A C
  //      C'   0];
  // where A'*A, or ATA, is given as a symmetric positive definite matrix and C
  // has full column-rank(?)
  //
  // [J] = lu_lagrange(ATA,C)
  //
  // Templates:
  //   T  should be a eigen matrix primitive type like int or double
  // Inputs:
  //   ATA   n by n square, symmetric, positive-definite system matrix, usually
  //     the quadratic coefficients corresponding to the original unknowns in a
  //     system
  //   C  n by m rectangular matrix corresponding the quadratic coefficients of
  //     the original unknowns times the lagrange multipliers enforcing linear
  //     equality constraints
  // Outputs:
  //   L  lower triangular matrix such that Q = L*U
  //   U  upper triangular matrix such that Q = L*U
  // Returns true on success, false on error
  //
  // Note: C should *not* have any empty columns. Typically C is the slice of
  // the linear constraints matrix Aeq concerning the unknown variables of a
  // quadratic optimization. Generally constraints may deal with unknowns as
  // well as knowns. Each linear constraint corresponds to a column of Aeq. As
  // long as each constraint concerns at least one unknown then the
  // corresponding column in C will have at least one non zero entry. If a
  // constraint concerns *no* unknowns, you should double check that this is a
  // valid constraint. How can you constrain known values to each other? This
  // is either a contradiction to the knowns' values or redundent. In either
  // case, it's not this functions responsiblilty to handle empty constraints
  // so you will get an error.
  //
  template <typename T>
  IGL_INLINE bool lu_lagrange(
    const Eigen::SparseMatrix<T> & ATA,
    const Eigen::SparseMatrix<T> & C,
    Eigen::SparseMatrix<T> & L,
    Eigen::SparseMatrix<T> & U);

}

#ifndef IGL_STATIC_LIBRARY
#  include "lu_lagrange.cpp"
#endif

#endif
