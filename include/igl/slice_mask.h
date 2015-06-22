// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SLICE_MASK_H
#define IGL_SLICE_MASK_H
#include "igl_inline.h"

#include <Eigen/Sparse>
#include <Eigen/Core>
namespace igl
{
  // Act like the matlab X(row_mask,col_mask) operator, where
  // row_mask, col_mask are non-negative integer indices.
  // 
  // Inputs:
  //   X  m by n matrix
  //   R  m list of row bools
  //   C  n list of column bools
  // Output:
  //   Y  #trues-in-R by #trues-in-C matrix
  //
  // See also: slice_mask
  
  template <typename DerivedX>
  IGL_INLINE void slice_mask(
    const Eigen::PlainObjectBase<DerivedX> & X,
    const Eigen::Array<bool,Eigen::Dynamic,1> & R,
    const Eigen::Array<bool,Eigen::Dynamic,1> & C,
    Eigen::PlainObjectBase<DerivedX> & Y);
  template <typename DerivedX>
  IGL_INLINE void slice_mask(
    const Eigen::PlainObjectBase<DerivedX> & X,
    const Eigen::Array<bool,Eigen::Dynamic,1> & R,
    const int dim,
    Eigen::PlainObjectBase<DerivedX> & Y);

  template <typename DerivedX>
  IGL_INLINE Eigen::PlainObjectBase<DerivedX> slice_mask(
    const Eigen::PlainObjectBase<DerivedX> & X,
    const Eigen::Array<bool,Eigen::Dynamic,1> & R,
    const Eigen::Array<bool,Eigen::Dynamic,1> & C);
  template <typename DerivedX>
  IGL_INLINE Eigen::PlainObjectBase<DerivedX> slice_mask(
    const Eigen::PlainObjectBase<DerivedX> & X,
    const Eigen::Array<bool,Eigen::Dynamic,1> & R,
    const int dim);
}


#ifndef IGL_STATIC_LIBRARY
#  include "slice_mask.cpp"
#endif

#endif
