#ifndef IGL_CROSS_H
#define IGL_CROSS_H
#include "igl_inline.h"
namespace igl
{
  // Computes out = cross(a,b)
  // Inputs:
  //   a  left 3d vector
  //   b  right 3d vector
  // Outputs:
  //   out  result 3d vector
  IGL_INLINE void cross(
    const double *a, 
    const double *b,
    double *out);
}

#ifdef IGL_HEADER_ONLY
#  include "cross.cpp"
#endif

#endif
