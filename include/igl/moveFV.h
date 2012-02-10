//
//  moveFV.h
//  Preview3D
//
//  Created by Olga Diamanti on 11/11/11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#ifndef IGL_MOVEFV_H
#define IGL_MOVEFV_H
#include "igl_inline.h"

#include <Eigen/Dense>
namespace igl 
{
  // moveFV 
  // Move a scalar field defined on faces to vertices by averaging
  //
  // Input:
  // V,F: mesh
  // S: scalar field defined on faces, Fx1
  // 
  // Output:
  // SV: scalar field defined on vertices
  template <typename T, typename I>
  IGL_INLINE void moveFV(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &V,
              const Eigen::Matrix<I, Eigen::Dynamic, Eigen::Dynamic> &F,
              const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &S,
              Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &SV);
}

#ifdef IGL_HEADER_ONLY
#  include "moveFV.cpp"
#endif

#endif
