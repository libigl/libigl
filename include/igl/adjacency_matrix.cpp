// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "adjacency_matrix.h"

#include "verbose.h"

// Bug in unsupported/Eigen/SparseExtra needs iostream first
#include <iostream>
#include <unsupported/Eigen/SparseExtra>

template <typename T>
IGL_INLINE void igl::adjacency_matrix(
  const Eigen::MatrixXi & F, 
  Eigen::SparseMatrix<T>& A)
{
  Eigen::DynamicSparseMatrix<T, Eigen::RowMajor> 
    dyn_A(F.maxCoeff()+1, F.maxCoeff()+1);
  switch(F.cols())
  {
    case 3:
      dyn_A.reserve(6*(F.maxCoeff()+1));
      break;
    case 4:
      dyn_A.reserve(26*(F.maxCoeff()+1));
      break;
  }

  // Loop over faces
  for(int i = 0;i<F.rows();i++)
  {
    // Loop over this face
    for(int j = 0;j<F.cols();j++)
    {
      // Get indices of edge: s --> d
      int s = F(i,j);
      int d = F(i,(j+1)%F.cols());
      dyn_A.coeffRef(s, d) = 1;
      dyn_A.coeffRef(d, s) = 1;
    }
  }

  A = Eigen::SparseMatrix<T>(dyn_A);
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template void igl::adjacency_matrix<int>(Eigen::Matrix<int, -1, -1, 0, -1, -1> const&, Eigen::SparseMatrix<int, 0, int>&);
#endif
