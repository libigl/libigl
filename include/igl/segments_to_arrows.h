// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Gustavo Segovia <gustavo.segovia@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SEGMENTS_TO_ARROWS_H
#define IGL_SEGMENTS_TO_ARROWS_H
#include "igl_inline.h"

#include <Eigen/Dense>

namespace igl
{
	
  // Takes segments from P1 to P2 and outputs segments that represent
  // an arrow pointing from P1 to P2.
  // Inputs:
  //   P1 Segment starting points
  //   P2 Segment end points
  // Outputs:
  //   A1 Arrow segment starting points
  //   A2 Arrow segment end points
  //
  // See also: adjacency_matrix
  IGL_INLINE void segments_to_arrows(const Eigen::MatrixXd& P1, const Eigen::MatrixXd& P2, Eigen::MatrixXd& A1, Eigen::MatrixXd& A2);
}

#ifndef IGL_STATIC_LIBRARY
#  include "segments_to_arrows.cpp"
#endif

#endif


