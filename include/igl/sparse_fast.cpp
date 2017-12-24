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
#include <unordered_map>
#include <map>
#include <utility>

// Hashing function
// namespace std {
//     template<> 
//     struct hash<pair<int, int> > {
//         size_t operator()(const pair<int, int>& p) const {
//             size_t seed = 0;
//             hash<int> h;
//             seed ^= h(p.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2); 
//             seed ^= h(p.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//             return seed;
//         }
//     };
// }


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

    // Build hash table for all nnz entries
    // TODO: this is slow and could be done in nlogn
    std::map<std::pair<int,int>,int> id;

    for (unsigned k=0; k<X.outerSize(); ++k)
    {
      unsigned outer_index = *(X.outerIndexPtr()+k);
      unsigned next_outer_index = (k+1 == X.outerSize()) ? X.nonZeros() : *(X.outerIndexPtr()+k+1); 
      
      for (unsigned l=outer_index; l<next_outer_index; ++l)
      {
        int col = k;
        int row = *(X.innerIndexPtr()+l);
        int value_index = l;

        std::pair<int,int> rc = std::make_pair(row,col);
        id[rc] = value_index;
      }
    }

    // Compute the indices
    data.resize(triplets.size());
    for (unsigned i=0; i<triplets.size(); ++i)
      data[i] = id[std::make_pair(triplets[i].row(),triplets[i].col())];
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
