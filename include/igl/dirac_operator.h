// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2018 Zhongshi Jiang <jiangzs@nyu.edu>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_DIRAC_H
#define IGL_DIRAC_H
#include <igl/igl_inline.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace igl 
{
  // Constructs the dirac operator with 4x4 matrix in place of quanternions
  // for a given mesh (V,F).
  //
  // Templates:
  //   DerivedV  derived type of eigen matrix for V (e.g. derived from
  //     MatrixXd)
  //   DerivedF  derived type of eigen matrix for F (e.g. derived from
  //     MatrixXi)
  //   Scalar  scalar type for eigen sparse matrix (e.g. double)
  // Inputs:
  //   V  #V by dim list of mesh vertex positions
  //   F  #F by simplex_size list of mesh faces (must be triangles)
  // Outputs: 
  //   D  4*#F by 4*#V dirac matrix
  //   DA  4*#V by 4*#F adjoint to dirac matrix
  //
  //
  // Note: This Dirac matrix uses a 4x4 matrix in place of quanternions
  //
  template <typename DerivedV, typename DerivedF, typename Scalar>
  IGL_INLINE void dirac_operator(
    const Eigen::MatrixBase<DerivedV> & V, 
    const Eigen::MatrixBase<DerivedF> & F,
    Eigen::SparseMatrix<Scalar>& D,
    Eigen::SparseMatrix<Scalar>& DA);
}

#ifndef IGL_STATIC_LIBRARY
#  include "dirac_operator.cpp"
#endif
#endif
