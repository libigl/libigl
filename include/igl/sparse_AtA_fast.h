// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2017 Daniele Panozzo <daniele.panozzo@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SPARSE_ATA_FAST_H
#define IGL_SPARSE_ATA_FAST_H
#include "igl_inline.h"
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>
namespace igl
{  
  struct sparse_AtA_fast_data
  {
    // Weights
    Eigen::VectorXd W;

    // Flatten composition rules
    std::vector<int> I_row;
    std::vector<int> I_col;
    std::vector<int> I_w;

    // For each entry of AtA, points to the beginning
    // of the composition rules
    std::vector<int> I_outer;
  };

  IGL_INLINE void sparse_AtA_fast_precompute(
    const Eigen::SparseMatrix<double>& A,
    Eigen::SparseMatrix<double>& AtA,
    sparse_AtA_fast_data& data);

  IGL_INLINE void sparse_AtA_fast(
    const Eigen::SparseMatrix<double>& A,
    Eigen::SparseMatrix<double>& AtA,
    const sparse_AtA_fast_data& data);
  
}

#ifndef IGL_STATIC_LIBRARY
#  include "sparse_AtA_fast.cpp"
#endif

#endif
