#ifndef IGL_UNIQUE_SIMPLICES_H
#define IGL_UNIQUE_SIMPLICES_H
#include "igl_inline.h"
#include <Eigen/Dense>
namespace igl
{
  // Find *combinatorially* unique simplices in F
  //
  // Inputs:
  //   F  #F by simplex-size list of simplices
  // Outputs:
  //   FF  #FF by simplex-size list of unique simplices in F
  IGL_INLINE void unique_simplices(
    const Eigen::MatrixXi & F,
    Eigen::MatrixXi & FF);
}

#ifdef IGL_HEADER_ONLY
#  include "unique_simplices.cpp"
#endif

#endif
