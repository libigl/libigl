#ifndef IGL_VIEW_AXIS_H
#define IGL_VIEW_AXIS_H 
#include "igl_inline.h"

namespace igl
{
  // Determines the view axis or depth axis of the current gl matrix
  // Outputs:
  //   x  pointer to x-coordinate in scene coordinates of the un-normalized
  //     viewing axis 
  //   y  pointer to y-coordinate in scene coordinates of the un-normalized
  //     viewing axis 
  //   z  pointer to z-coordinate in scene coordinates of the un-normalized
  //     viewing axis
  //   mv pointer to modelview matrix
  //
  // Note: View axis is returned *UN-normalized*
  IGL_INLINE void view_axis(double * x, double * y, double * z);
  IGL_INLINE void view_axis(const double * mv, double * x, double * y, double * z);
};


#ifdef IGL_HEADER_ONLY
#  include "view_axis.cpp"
#endif

#endif

