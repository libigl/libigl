// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "diag.h"

#include "verbose.h"

// Bug in unsupported/Eigen/SparseExtra needs iostream first
#include <iostream>
#include <unsupported/Eigen/SparseExtra>

template <typename T>
IGL_INLINE void igl::diag(
  const Eigen::SparseMatrix<T>& X,
  Eigen::SparseVector<T>& V)
{
  assert(false && "Just call X.diagonal().sparseView() directly");
  V = X.diagonal().sparseView();
}

template <typename T,typename DerivedV>
IGL_INLINE void igl::diag(
  const Eigen::SparseMatrix<T>& X,
  Eigen::MatrixBase<DerivedV> & V)
{
  assert(false && "Just call X.diagonal() directly");
  V = X.diagonal();
}

template <typename T>
IGL_INLINE void igl::diag(
  const Eigen::SparseVector<T>& V,
  Eigen::SparseMatrix<T>& X)
{
  // clear and resize output
  std::vector<Eigen::Triplet<double> > Xijv;
  const int n = V.size();
  const int nnz = V.nonZeros();
  Xijv.reserve(nnz);
  // loop over non-zeros
  for(typename Eigen::SparseVector<T>::InnerIterator it(V); it; ++it)
  {
    Xijv.emplace_back(it.index(),it.index(),it.value());
  }
  X.resize(n,n);
  X.setFromTriplets(Xijv.begin(),Xijv.end());
}

template <typename T, typename DerivedV>
IGL_INLINE void igl::diag(
  const Eigen::MatrixBase<DerivedV> & V,
  Eigen::SparseMatrix<T>& X)
{
  assert(V.rows() == 1 || V.cols() == 1);
  // clear and resize output
  std::vector<Eigen::Triplet<double> > Xijv;
  const int n = V.size();
  Xijv.reserve(n);
  // loop over non-zeros
  for(int i = 0;i<n;i++) { Xijv.emplace_back(i,i,V(i)); }
  X.resize(n,n);
  X.setFromTriplets(Xijv.begin(),Xijv.end());
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::diag<double, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::SparseMatrix<double, 0, int> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::diag<double>(Eigen::SparseMatrix<double, 0, int> const&, Eigen::SparseVector<double, 0, int>&);
template void igl::diag<double, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::SparseMatrix<double, 0, int>&);
template void igl::diag<double>(Eigen::SparseVector<double, 0, int> const&, Eigen::SparseMatrix<double, 0, int>&);
#endif
