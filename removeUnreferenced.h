//
//  removeUnreferenced.h
//  Preview3D
//
//  Created by Daniele Panozzo on 17/11/11.

#ifndef IGL_REMOVEUNREFERENCED_H
#define IGL_REMOVEUNREFERENCED_H
#include "igl_inline.h"

#include <Eigen/Core>
namespace igl 
{
  // [ NV, NF ] = removeUnreferenced( V,F,epsilon )
  // Remove unreferenced vertices from V, updating F accordingly
  //
  // Input:
  // V,F: mesh description
  //
  // Output:
  // NV, NF: new mesh without unreferenced vertices
  
  template <typename T>
  IGL_INLINE void removeUnreferenced(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &V, const Eigen::MatrixXi &F, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &NV, Eigen::MatrixXi &NF, Eigen::VectorXi &I);
  
}

#ifdef IGL_HEADER_ONLY
#  include "removeUnreferenced.cpp"
#endif

#endif
