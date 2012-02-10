//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

#ifndef IGL_IS_BORDER_VERTEX_H
#define IGL_IS_BORDER_VERTEX_H
#include "igl_inline.h"

#include <Eigen/Core>
#include <vector>

namespace igl 
{
  template<typename T, typename S>
  IGL_INLINE std::vector<bool> is_border_vertex(const T& V,
                                                const Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>& F);
}

#ifdef IGL_HEADER_ONLY
#  include "is_border_vertex.cpp"
#endif

#endif
