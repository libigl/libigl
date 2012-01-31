//
//  removeDuplicates.h
//  Preview3D
//
//  Created by Olga Diamanti on 17/11/11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#ifndef IGL_REMOVEDUPLICATES_H
#define IGL_REMOVEDUPLICATES_H
#include "igl_inline.h"

#include <Eigen/Core>
namespace igl 
{
  // [ NV, NF ] = removeDuplicates( V,F,epsilon )
  // Merge the duplicate vertices from V, fixing the topology accordingly
  //
  // Input:
  // V,F: mesh description
  // epsilon: minimal distance to consider two vertices identical
  //
  // Output:
  // NV, NF: new mesh without duplicate vertices
  
  template <typename T>
  IGL_INLINE void removeDuplicates(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &V, const Eigen::MatrixXi &F, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &NV, Eigen::MatrixXi &NF, Eigen::VectorXi &I, const double epsilon = 2.2204e-15);
  
}

#ifdef IGL_HEADER_ONLY
#  include "removeDuplicates.cpp"
#endif

#endif
