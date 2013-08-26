#ifndef IGL_RIGHT_AXIS_H
#define IGL_RIGHT_AXIS_H 
#include "igl_inline.h"
namespace igl
{
  // Determines the right axis or depth axis of the current gl matrix
  // Outputs:
  //   x  pointer to x-coordinate in scene coordinates of the un-normalized
  //     right axis 
  //   y  pointer to y-coordinate in scene coordinates of the un-normalized
  //     right axis 
  //   z  pointer to z-coordinate in scene coordinates of the un-normalized
  //     right axis
  //   mv pointer to modelview matrix
  //
  // Note: Right axis is returned *UN-normalized*
  IGL_INLINE void right_axis(double * x, double * y, double * z);
  IGL_INLINE void right_axis(const double * mv, double * x, double * y, double * z);
};

#ifdef IGL_HEADER_ONLY
#  include "right_axis.cpp"
#endif
#endif
