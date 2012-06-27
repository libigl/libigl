#ifndef IGL_DRAW_POINT_H
#define IGL_DRAW_POINT_H
#include "igl_inline.h"
namespace igl
{
  //double POINT_COLOR[3] = {239./255.,213./255.,46./255.};
  // Draw a nice looking, colored dot at a given point in 3d.
  //
  // Note: expects that GL_CURRENT_COLOR is set with the desired foreground color
  // 
  // Inputs:
  //   x  x-coordinate of point, modelview coordinates
  //   y  y-coordinate of point, modelview coordinates
  //   z  z-coordinate of point, modelview coordinates
  //   requested_r  outer-most radius of dot {7}, measured in screen space pixels
  //   selected  fills inner circle with black {false}
  // Asserts that requested_r does not exceed 0.5*GL_POINT_SIZE_MAX
  IGL_INLINE void draw_point(
    const double x,
    const double y,
    const double z,
    const double requested_r = 7,
    const bool selected = false);
}

#ifdef IGL_HEADER_ONLY
#  include "draw_point.cpp"
#endif

#endif
