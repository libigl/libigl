#ifndef IGL_DOT_H
#define IGL_DOT_H
#include "igl_inline.h"
namespace igl
{
  // Computes out = dot(a,b)
  // Inputs:
  //   a  left 3d vector
  //   b  right 3d vector
  // Returns scalar dot product
  IGL_INLINE double dot(
    const double *a, 
    const double *b);
}

#ifdef IGL_HEADER_ONLY
#  include "dot.cpp"
#endif

#endif
