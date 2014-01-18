// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_HARMONIC_H
#define IGL_HARMONIC_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl 
{
  // Compute k-harmonic weight functions "coordinates".
  //
  //
  // Inputs:
  //   V  #V by dim vertex positions
  //   F  #F by simplex-size list of element indices
  //   b  #b boundary indices into V
  //   bc #b by #W list of boundary values
  //   k  power of harmonic operation (1: harmonic, 2: biharmonic, etc)
  // Outputs:
  //   W  #V by #W list of weights
  //
  IGL_INLINE bool harmonic(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::VectorXi & b,
    const Eigen::MatrixXd & bc,
    const int k,
    Eigen::MatrixXd & W);
};
#ifdef IGL_HEADER_ONLY
#include "harmonic.cpp"
#endif
#endif
