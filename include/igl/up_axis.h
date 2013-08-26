#ifndef IGL_UP_AXIS_H
#define IGL_UP_AXIS_H 
#include "igl_inline.h"
namespace igl
{
// Determines the up axis or depth axis of the current gl matrix
// Outputs:
//   x  pointer to x-coordinate in scene coordinates of the un-normalized
//     up axis 
//   y  pointer to y-coordinate in scene coordinates of the un-normalized
//     up axis 
//   z  pointer to z-coordinate in scene coordinates of the un-normalized
//     up axis
  //   mv pointer to modelview matrix
//
// Note: Up axis is returned *UN-normalized*
IGL_INLINE void up_axis(double * x, double * y, double * z);
IGL_INLINE void up_axis(const double * mv, double * x, double * y, double * z);
};

#ifdef IGL_HEADER_ONLY
#  include "up_axis.cpp"
#endif
#endif


