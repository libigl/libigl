// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "slice_cached.h"

#include <iostream>
#include <vector>
#include <utility>
#include <igl/slice.h>

template <typename TX, typename TY>
IGL_INLINE void igl::slice_cached_precompute(
  const Eigen::SparseMatrix<TX>& X,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & C,
  Eigen::SparseMatrix<TY>& Y,
  Eigen::VectorXi& data)
{
  // Create a sparse matrix whose entries are the ids
  Eigen::SparseMatrix<unsigned> TS = X.template cast<unsigned>();

  TS.makeCompressed();
  for (unsigned i=0;i<TS.nonZeros();++i)
    *(TS.valuePtr() + i) = i;

  Eigen::SparseMatrix<unsigned> TS_sliced;
  igl::slice(TS,R,C,TS_sliced);
  Y = TS_sliced.cast<TY>();

  data.resize(TS_sliced.nonZeros());
  for (unsigned i=0;i<data.size();++i)
  {
    data[i] = *(TS_sliced.valuePtr() + i);
    *(Y.valuePtr() + i) = *(X.valuePtr() + data[i]);
  }
}

template <typename TX, typename TY>
IGL_INLINE void igl::slice_cached(
  const Eigen::SparseMatrix<TX>& X,
  Eigen::SparseMatrix<TY>& Y,
  const Eigen::VectorXi& data)
{
  for (unsigned i=0; i<data.size(); ++i)
    *(Y.valuePtr() + i) = *(X.valuePtr() + data[i]);
}

#ifdef IGL_STATIC_LIBRARY
template void igl::slice_cached_precompute<double, double>(Eigen::SparseMatrix<double, 0, int> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::SparseMatrix<double, 0, int>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&);
template void igl::slice_cached<double, double>(Eigen::SparseMatrix<double, 0, int> const&, Eigen::SparseMatrix<double, 0, int>&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&);
#endif
