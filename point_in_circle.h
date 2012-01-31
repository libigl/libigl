#ifndef IGL_POINT_IN_CIRCLE_H
#define IGL_POINT_IN_CIRCLE_H
#include "igl_inline.h"

namespace igl
{
  // Determine if 2d point is in a circle
  // Inputs:
  //   qx  x-coordinate of query point
  //   qy  y-coordinate of query point
  //   cx  x-coordinate of circle center
  //   cy  y-coordinate of circle center
  //   r  radius of circle
  // Returns true if query point is in circle, false otherwise
  IGL_INLINE bool point_in_circle(
    const double qx, 
    const double qy,
    const double cx, 
    const double cy,
    const double r);
}

#ifdef IGL_HEADER_ONLY
#  include "point_in_circle.cpp"
#endif

#endif
