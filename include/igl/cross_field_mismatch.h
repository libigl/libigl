// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>, Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_CROSS_FIELD_MISMATCH_H
#define IGL_CROSS_FIELD_MISMATCH_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Calculates the mismatch (integer), at each face edge, of a cross field defined on the mesh faces.
  // The integer mismatch is a multiple of pi/2 that transforms the cross on one side of the edge to
  // the cross on the other side. It represents the deviation from a Lie connection across the edge.

  // Inputs:
  //   V         #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F         #F by 3 eigen Matrix of face (quad) indices
  //   PD1       #F by 3 eigen Matrix of the first per face cross field vector
  //   PD2       #F by 3 eigen Matrix of the second per face cross field vector
  //   isCombed  boolean, specifying whether the field is combed (i.e. matching has been precomputed.
  //             If not, the field is combed first.
  // Output:
  //   mismatch  #F by 3 eigen Matrix containing the integer mismatch of the cross field
  //             across all face edges
  //

  template <typename DerivedV, typename DerivedF, typename DerivedM>
  IGL_INLINE void cross_field_mismatch(const Eigen::MatrixBase<DerivedV> &V,
                                       const Eigen::MatrixBase<DerivedF> &F,
                                       const Eigen::MatrixBase<DerivedV> &PD1,
                                       const Eigen::MatrixBase<DerivedV> &PD2,
                                       const bool isCombed,
                                       Eigen::PlainObjectBase<DerivedM> &mismatch);
}
#ifndef IGL_STATIC_LIBRARY
#include "cross_field_mismatch.cpp"
#endif

#endif
