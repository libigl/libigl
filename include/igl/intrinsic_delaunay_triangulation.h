// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2018 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_INTRINSIC_DELAUNAY_TRIANGULATION_H
#define IGL_INTRINSIC_DELAUNAY_TRIANGULATION_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // INTRINSIC_DELAUNAY_TRIANGULATION Flip edges _intrinsically_ until all are
  // "intrinsic Delaunay". See "An algorithm for the construction of intrinsic
  // delaunay triangulations with applications to digital geometry processing"
  // [Fisher et al. 2007].
  //
  // Inputs:
  //   l_in  #F_in by 3 list of edge lengths (see edge_lengths)
  //   F_in  #F_in by 3 list of face indices into some unspecified vertex list V
  // Outputs:
  //   l  #F by 3 list of edge lengths
  //   F  #F by 3 list of new face indices
  //
  // See also: is_intrinsic_delaunay
  template <
    typename Derivedl_in,
    typename DerivedF_in,
    typename Derivedl,
    typename DerivedF>
  IGL_INLINE void intrinsic_delaunay_triangulation(
    const Eigen::MatrixBase<Derivedl_in> & l_in,
    const Eigen::MatrixBase<DerivedF_in> & F_in,
    Eigen::PlainObjectBase<Derivedl> & l,
    Eigen::PlainObjectBase<DerivedF> & F);
}

#ifndef IGL_STATIC_LIBRARY
#  include "intrinsic_delaunay_triangulation.cpp"
#endif

#endif
