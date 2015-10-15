// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_FIELD_LOCAL_GLOBAL_CONVERSIONS
#define IGL_FIELD_LOCAL_GLOBAL_CONVERSIONS
#include "igl_inline.h"

#include <Eigen/Core>
#include <vector>

namespace igl {
  // Converts a face-based polyvector field consisting of n vectors per face
  // from its 3D coordinates to its local 2D representation (with respect to the
  // local bases of each triangle)

  // Inputs:
  //   B1               #F by 3 list of the first basis vector of each triangle
  //   B2               #F by 3 list of the second basis vector of each triangle
  //   global           #F by 3n list of the 3D coordinates of the per-face vectors
  //                    (stacked horizontally for each triangle)
  // Output:
  //   local            #F by 2n list of the 2D representation of the per-face vectors
  //                    (stacked horizontally for each triangle)
  //
template <typename DerivedG, typename DerivedL, typename DerivedB>
IGL_INLINE void global2local(
const Eigen::PlainObjectBase<DerivedB>& B1,
const Eigen::PlainObjectBase<DerivedB>& B2,
const Eigen::PlainObjectBase<DerivedG>& global,
Eigen::PlainObjectBase<DerivedL>& local);

// Converts a face-based polyvector field consisting of n vectors per face
// from its local 2D representation (with respect to the local bases of each
// triangle) to its 3D coordinates

// Inputs:
//   B1               #F by 3 list of the first basis vector of each triangle
//   B2               #F by 3 list of the second basis vector of each triangle
//   local            #F by 2n list of the 2D representation of the per-face vectors
//                    (stacked horizontally for each triangle)
// Output:
//   global           #F by 3n list of the 3D coordinates of the per-face vectors
//                    (stacked horizontally for each triangle)
//
template <typename DerivedG, typename DerivedL, typename DerivedB>
IGL_INLINE void local2global(
const Eigen::PlainObjectBase<DerivedB>& B1,
const Eigen::PlainObjectBase<DerivedB>& B2,
const Eigen::PlainObjectBase<DerivedL>& local,
Eigen::PlainObjectBase<DerivedG>& global);

};


#ifndef IGL_STATIC_LIBRARY
#include "field_local_global_conversions.cpp"
#endif


#endif /* defined(IGL_FIELD_LOCAL_GLOBAL_CONVERSIONS) */
