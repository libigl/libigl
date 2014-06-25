// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
//
//  remove_unreferenced.h
//  Preview3D
//
//  Created by Daniele Panozzo on 17/11/11.

#ifndef IGL_REMOVE_UNREFERENCED_H
#define IGL_REMOVE_UNREFERENCED_H
#include "igl_inline.h"

#include <Eigen/Core>
namespace igl 
{
  // [ NV, NF ] = remove_unreferenced( V,F)
  // Remove unreferenced vertices from V, updating F accordingly
  //
  // Input:
  // V,F: mesh description
  //
  // Output:
  // NV, NF: new mesh without unreferenced vertices
  //
  template <typename Scalar, typename Index>
  IGL_INLINE void remove_unreferenced(
    const Eigen::PlainObjectBase<Scalar> &V,
    const Eigen::PlainObjectBase<Index> &F,
    Eigen::PlainObjectBase<Scalar> &NV,
    Eigen::PlainObjectBase<Index> &NF,
    Eigen::PlainObjectBase<Index> &I);
}

#ifndef IGL_STATIC_LIBRARY
#  include "remove_unreferenced.cpp"
#endif

#endif
