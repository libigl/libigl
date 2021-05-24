// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_WINDING_NUMBER_H
#define IGL_WINDING_NUMBER_H
#include "igl_inline.h"
#include <Eigen/Core>

// Minimum number of iterms per openmp thread
#ifndef IGL_WINDING_NUMBER_OMP_MIN_VALUE
#  define IGL_WINDING_NUMBER_OMP_MIN_VALUE 1000
#endif
namespace igl
{
  // WINDING_NUMBER  Computes the generalized winding number at each
  // dim-dimensional query point in O with respect to the oriented
  // one-codimensional mesh (V,F). This is equivalent to summing the subtended
  // signed angles/solid angles of each element in (V,F). See, "Robust
  // Inside-Outside Segmentation using Generalized Winding Numbers" [Jacobson et
  // al. 2013].
  //
  //
  // Inputs:
  //   V  #V by dim list of mesh vertex positions
  //   F  #F by dim list of mesh facets as indices into rows of V. If dim==2,
  //     then (V,F) describes a set of edges in the plane. If dim==3, then (V,F)
  //     describes a triangle mesh/soup.
  //   O  #O by dim list of query points
  // Output:
  //   W  #O by 1 list of winding numbers
  //
  // See also: igl::fast_winding_number
  //
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedO,
    typename DerivedW>
  IGL_INLINE void winding_number(
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    const Eigen::MatrixBase<DerivedO> & O,
    Eigen::PlainObjectBase<DerivedW> & W);
  // Compute winding number of a single point
  //
  // Inputs:
  //  V  n by dim list of vertex positions
  //  F  #F by dim list of triangle indices, minimum index is 0
  //  p  single origin position
  // Outputs:
  //  w  winding number of this point
  //
  template <
    typename DerivedV,
    typename DerivedF,
    typename Derivedp>
  IGL_INLINE typename DerivedV::Scalar winding_number(
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    const Eigen::MatrixBase<Derivedp> & p);
}

#ifndef IGL_STATIC_LIBRARY
#  include "winding_number.cpp"
#endif

#endif
