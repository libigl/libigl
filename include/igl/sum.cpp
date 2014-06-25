// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "sum.h"

template <typename T>
IGL_INLINE void igl::sum(
  const Eigen::SparseMatrix<T>& X, 
  const int dim,
  Eigen::SparseVector<T>& S)
{
  // dim must be 2 or 1
  assert(dim == 1 || dim == 2);
  // Get size of input
  int m = X.rows();
  int n = X.cols();
  // resize output
  if(dim==1)
  {
    S = Eigen::SparseVector<T>(n);
  }else
  {
    S = Eigen::SparseVector<T>(m);
  }

  // Iterate over outside
  for(int k=0; k<X.outerSize(); ++k)
  {
    // Iterate over inside
    for(typename Eigen::SparseMatrix<T>::InnerIterator it (X,k); it; ++it)
    {
      if(dim == 1)
      {
        S.coeffRef(it.col()) += it.value();
      }else
      {
        S.coeffRef(it.row()) += it.value();
      }
    }
  }

}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
#endif
