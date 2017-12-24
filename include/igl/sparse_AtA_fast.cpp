// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "sparse_AtA_fast.h"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <map>
#include <utility>

IGL_INLINE void igl::sparse_AtA_fast_precompute(
    const Eigen::SparseMatrix<double>& A,
    Eigen::SparseMatrix<double>& AtA,
    igl::sparse_AtA_fast_data& data)
{
  // 1 Compute At
  std::vector<std::vector<int> > Col_RowPtr;
  std::vector<std::vector<int> > Col_IndexPtr;

  Col_RowPtr.resize(A.cols());
  Col_IndexPtr.resize(A.cols());

  for (unsigned k=0; k<A.outerSize(); ++k)
  {
    unsigned outer_index = *(A.outerIndexPtr()+k);
    unsigned next_outer_index = (k+1 == A.outerSize()) ? A.nonZeros() : *(A.outerIndexPtr()+k+1); 
    
    for (unsigned l=outer_index; l<next_outer_index; ++l)
    {
      int col = k;
      int row = *(A.innerIndexPtr()+l);
      int value_index = l;
      assert(col < A.cols());
      assert(col >= 0);
      assert(row < A.rows());
      assert(row >= 0);
      assert(value_index >= 0);
      assert(value_index < A.nonZeros());

      Col_RowPtr.at(col).push_back(row);
      Col_IndexPtr.at(col).push_back(value_index);
    }
  }

  Eigen::SparseMatrix<double> At = A.transpose();
  At.makeCompressed();
  AtA = At * A;
  AtA.makeCompressed();

  assert(AtA.isCompressed());

  // If weights are not provided, use 1
  if (data.W.size() == 0)
    data.W = Eigen::VectorXd::Ones(A.rows());
  assert(data.W.size() == A.rows());

  // 2 Construct the rules
  for (unsigned k=0; k<AtA.outerSize(); ++k)
  {
    unsigned outer_index = *(AtA.outerIndexPtr()+k);
    unsigned next_outer_index = (k+1 == AtA.outerSize()) ? AtA.nonZeros() : *(AtA.outerIndexPtr()+k+1); 
    
    for (unsigned l=outer_index; l<next_outer_index; ++l)
    {
      int col = k;
      int row = *(AtA.innerIndexPtr()+l);
      int value_index = l;
      assert(col < AtA.cols());
      assert(col >= 0);
      assert(row < AtA.rows());
      assert(row >= 0);
      assert(value_index >= 0);
      assert(value_index < AtA.nonZeros());

      data.I_outer.push_back(data.I_row.size());

      // Find corresponding indices
      for (unsigned i=0;i<Col_RowPtr.at(row).size();++i)
      {
        for (unsigned j=0;j<Col_RowPtr.at(col).size();++j)
        {
          if (Col_RowPtr.at(row)[i] == Col_RowPtr.at(col)[j])
          {
            data.I_row.push_back(Col_IndexPtr[row][i]);
            data.I_col.push_back(Col_IndexPtr[col][j]);
            data.I_w.push_back(data.W[Col_RowPtr[col][j]]);
          }
        }
      }
    }
  }
  data.I_outer.push_back(data.I_row.size()); // makes it more efficient to iterate later on
}

IGL_INLINE void igl::sparse_AtA_fast(
    const Eigen::SparseMatrix<double>& A,
    Eigen::SparseMatrix<double>& AtA,
    const igl::sparse_AtA_fast_data& data)
{
  for (unsigned i=0; i<data.I_outer.size()-1; ++i)
  {
    *(AtA.valuePtr() + i) = 0;
    for (unsigned j=data.I_outer[i]; j<data.I_outer[i+1]; ++j)
      *(AtA.valuePtr() + i) += *(A.valuePtr() + data.I_row[j]) * data.I_w[j] * *(A.valuePtr() + data.I_col[j]);
  }
}


#ifdef IGL_STATIC_LIBRARY
#endif
