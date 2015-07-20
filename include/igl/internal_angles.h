// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_INTERNAL_ANGLES_H
#define IGL_INTERNAL_ANGLES_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Compute internal angles for a triangle mesh
  //
  // Inputs:
  //   V  #V by dim eigen Matrix of mesh vertex nD positions
  //   F  #F by poly-size eigen Matrix of face (triangle) indices
  // Output:
  //   K  #F by poly-size eigen Matrix of internal angles
  //     for triangles, columns correspond to edges [1,2],[2,0],[0,1]
  //
  // Known Issues:
  //   if poly-size ≠ 3 then dim must equal 3.
  template <typename DerivedV, typename DerivedF, typename DerivedK>
  IGL_INLINE void internal_angles(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedK> & K);
  // Inputs:
  //   L  #F by 3 list of edge lengths
  template <typename DerivedL, typename DerivedK>
  IGL_INLINE void internal_angles(
    const Eigen::PlainObjectBase<DerivedL>& L,
    Eigen::PlainObjectBase<DerivedK> & K);
}

#ifndef IGL_STATIC_LIBRARY
#  include "internal_angles.cpp"
#endif

#endif


