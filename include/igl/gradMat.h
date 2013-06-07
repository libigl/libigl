//
//  gradMat.h
//
//  Created by Christian Schüller on 07/03/13.
//  Copyright (c) 2013 ETHZ IGL. All rights reserved.
//

#ifndef IGL_GRAD_MAT_H
#define IGL_GRAD_MAT_H
#include "igl_inline.h"

#include <Eigen/Core>

namespace igl {
  // GRAD
  // G = grad(V,F)
  //
  // Compute the numerical gradient operator
  //
  // Inputs:
  //   V  #vertices by 3 list of mesh vertex positions
  //   F  #faces by 3 list of mesh face indices
  // Outputs:
  //   G  #faces*dim by #V Gradient operator
  //
  
  // Gradient of a scalar function defined on piecewise linear elements (mesh)
  // is constant on each triangle i,j,k:
  // grad(Xijk) = (Xj-Xi) * (Vi - Vk)^R90 / 2A + (Xk-Xi) * (Vj - Vi)^R90 / 2A
  // where Xi is the scalar value at vertex i, Vi is the 3D position of vertex
  // i, and A is the area of triangle (i,j,k). ^R90 represent a rotation of 
  // 90 degrees
  //
  template <typename T, typename S>
  IGL_INLINE void gradMat(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &V,
    const Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic> &F,
  Eigen::SparseMatrix<T> &G);
}

#ifdef IGL_HEADER_ONLY
#  include "gradMat.cpp"
#endif

#endif
