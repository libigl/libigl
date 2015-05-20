// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>, Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_ROTATION_MATRIX_FROM_DIRECTIONS
#define IGL_ROTATION_MATRIX_FROM_DIRECTIONS
#include "igl_inline.h"

#include <Eigen/Core>

namespace igl {
  /// Given 2 vectors centered on origin calculate the rotation matrix from first to the second

  // Inputs:
  //   v0, v1         the two #3 by 1 vectors
  //   normalized     boolean, if false, then the vectors are normalized prior to the calculation
  // Output:
  //                  3 by 3 rotation matrix that takes v0 to v1
  //
  template <typename Scalar>
  IGL_INLINE Eigen::Matrix<Scalar, 3, 3> rotation_matrix_from_directions(const Eigen::Matrix<Scalar, 3, 1> v0,
                                                                     const Eigen::Matrix<Scalar, 3, 1> v1,
                                                                     const bool normalized=true);
}


#ifndef IGL_STATIC_LIBRARY
#include "rotation_matrix_from_directions.cpp"
#endif


#endif /* defined(IGL_ROTATION_MATRIX_FROM_DIRECTIONS) */
