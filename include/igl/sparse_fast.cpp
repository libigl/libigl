// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "sparse_fast.h"

#include <iostream>
#include <vector>
#include <array>
#include <unordered_map>
#include <map>
#include <utility>

IGL_INLINE void igl::sparse_fast_precompute(
  const Eigen::VectorXi & I,
  const Eigen::VectorXi & J,
  Eigen::SparseMatrix<double>& X,
  Eigen::VectorXi& data)
{
  // Generates the triplets
  std::vector<Eigen::Triplet<double> > t(I.size());
  for (unsigned i = 0; i<I.size(); ++i)
    t[i] = Eigen::Triplet<double>(I[i],J[i],1);

  // Call the triplets version
  sparse_fast_precompute(t,X,data);
}

IGL_INLINE void igl::sparse_fast_precompute(
  const std::vector<Eigen::Triplet<double> >& triplets,
  Eigen::SparseMatrix<double>& X,
  Eigen::VectorXi& data)
{
    // Construct an empty sparse matrix
    X.setFromTriplets(triplets.begin(),triplets.end());
    X.makeCompressed();

    std::vector<std::array<int,3> > T(triplets.size());
    for (unsigned i=0; i<triplets.size(); ++i)
    {
      T[i][0] = triplets[i].col();
      T[i][1] = triplets[i].row();
      T[i][2] = i;
    }

    std::sort(T.begin(), T.end());

    data.resize(triplets.size());

    int t = 0;

    for (unsigned k=0; k<X.outerSize(); ++k)
    {
      unsigned outer_index = *(X.outerIndexPtr()+k);
      unsigned next_outer_index = (k+1 == X.outerSize()) ? X.nonZeros() : *(X.outerIndexPtr()+k+1); 
      
      for (unsigned l=outer_index; l<next_outer_index; ++l)
      {
        int col = k;
        int row = *(X.innerIndexPtr()+l);
        int value_index = l;
        assert(col < X.cols());
        assert(col >= 0);
        assert(row < X.rows());
        assert(row >= 0);
        assert(value_index >= 0);
        assert(value_index < X.nonZeros());

        std::pair<int,int> p_m = std::make_pair(row,col);

        while (t<T.size() && (p_m == std::make_pair(T[t][1],T[t][0])))
          data[T[t++][2]] = value_index;
      }
    }
    assert(t==T.size());

}
  
IGL_INLINE void igl::sparse_fast(
  const std::vector<Eigen::Triplet<double> >& triplets,
  Eigen::SparseMatrix<double>& X,
  const Eigen::VectorXi& data)
{
  assert(triplets.size() == data.size());

  // Clear it first
  for (unsigned i = 0; i<data.size(); ++i)
    *(X.valuePtr() + data[i]) = 0;
 
  // Then sum them up
  for (unsigned i = 0; i<triplets.size(); ++i)
    *(X.valuePtr() + data[i]) += triplets[i].value();
}

IGL_INLINE void igl::sparse_fast(
  const Eigen::VectorXd & V,
  Eigen::SparseMatrix<double>& X,
  const Eigen::VectorXi& data)
{
  assert(V.size() == data.size());

  // Clear it first
  for (unsigned i = 0; i<data.size(); ++i)
    *(X.valuePtr() + data[i]) = 0;
 
  // Then sum them up
  for (unsigned i = 0; i<V.size(); ++i)
    *(X.valuePtr() + data[i]) += V[i];
}


#ifdef IGL_STATIC_LIBRARY
#endif
