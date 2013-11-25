// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_EDGETOPOLOGY_H
#define IGL_EDGETOPOLOGY_H

#include "igl_inline.h"

#include <Eigen/Core>
#include <vector>

namespace igl 
{
  // Initialize Edges and their topological relations
  
  // Output:
  // EV  : #Ex2, Stores the edge description as pair of indices to vertices
  // FE : #Fx3, Stores the Triangle-Edge relation
  // EF : #Ex2: Stores the Edge-Triangle relation
  IGL_INLINE void edgetopology(
    const Eigen::MatrixXd& V, 
    const Eigen::MatrixXi& F, 
    Eigen::MatrixXi& EV, 
    Eigen::MatrixXi& FE, 
    Eigen::MatrixXi& EF);
}

#ifdef IGL_HEADER_ONLY
#  include "edgetopology.cpp"
#endif

#endif
