//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

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
  // EF : #Ex2: Stores the Edge-Triangle relation (unsorted)
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
