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
  // WINDING_NUMBER Compute the sum of solid angles of a triangle/tetrahedron
  // described by points (vectors) V
  //
  // Templates:
  //   dim  dimension of input
  // Inputs:
  //  V  n by 3 list of vertex positions
  //  F  #F by 3 list of triangle indices, minimum index is 0
  //  O  no by 3 list of origin positions
  // Outputs:
  //  S  no by 1 list of winding numbers
  //
  // 3d
  IGL_INLINE void winding_number(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & O,
    Eigen::VectorXd & W);
  // Inputs:
  //   V  pointer to array containing #V by 3 vertex positions along rows,
  //     given in column major order
  //   n  number of mesh vertices
  //   F  pointer to array containing #F by 3 face indices along rows,
  //     given in column major order
  //   m  number of faces
  //   O  pointer to array containing #O by 3 query positions along rows,
  //     given in column major order
  //   no  number of origins
  // Outputs:
  //   S  no by 1 list of winding numbers
  template <typename Scalar, typename DerivedF>
  IGL_INLINE void winding_number_3(
    const Scalar * V,
    const int n,
    const DerivedF * F,
    const int m,
    const Scalar * O,
    const int no,
    Scalar * S);
  //// Only one evaluation origin
  //template <typename DerivedF>
  //IGL_INLINE void winding_number_3(
  //  const double * V,
  //  const int n,
  //  const DerivedF * F,
  //  const int m,
  //  const double * O,
  //  double * S);
  // 2d
  template <typename DerivedF>
  IGL_INLINE void winding_number_2(
    const double * V,
    const int n,
    const DerivedF * F,
    const int m,
    const double * O,
    const int no,
    double * S);
}

#ifndef IGL_STATIC_LIBRARY
#  include "winding_number.cpp"
#endif

#endif
