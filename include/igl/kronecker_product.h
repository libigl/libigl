// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_KRONECKERPRODUCT_H
#define IGL_KRONECKERPRODUCT_H
#include "igl_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace igl
{
  // Computes the Kronecker product between sparse matrices A and B.
  //
  // Inputs:
  //   A  #M by #N sparse matrix
  //   B  #P by #Q sparse matrix
  // Outputs:
  //      #M*#P by #N*#Q sparse matrix
  //
  template <typename Scalar>
  IGL_INLINE Eigen::SparseMatrix<Scalar> kronecker_product(
    const Eigen::SparseMatrix<Scalar> & A,
    const Eigen::SparseMatrix<Scalar> & B);
}

#ifdef IGL_HEADER_ONLY
#  include "kronecker_product.cpp"
#endif

#endif
