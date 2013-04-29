#ifndef IGL_MVC_H
#define IGL_MVC_H

#include "igl_inline.h"
#include <Eigen/Dense>

namespace igl 
{
  //   MVC - MEAN VALUE COORDINATES
  //  
  //   mvc(V,C,W)
  //  
  //   Inputs:
  //    V  #V x dim list of vertex positions (dim = 2 or dim = 3)
  //    C  #C x dim list of polygon vertex positions in counter-clockwise order
  //      (dim = 2 or dim = 3)
  //  
  //   Outputs:
  //    W  weights, #V by #C matrix of weights
  //  
  //  Known Bugs: implementation is listed as "Broken"
  IGL_INLINE void mvc(
    const Eigen::MatrixXd &V, 
    const Eigen::MatrixXd &C, 
    Eigen::MatrixXd &W);
  
}

#ifdef IGL_HEADER_ONLY
#  include "mvc.cpp"
#endif

#endif
