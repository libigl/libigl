#ifndef IGL_RGB_TO_HSV_H
#define IGL_RGB_TO_HSV_H
#include "igl_inline.h"
namespace igl
{
  // Convert RGB to HSV
  //
  // Inputs:
  //   r  red value ([0,1]) 
  //   g  green value ([0,1])
  //   b  blue value ([0,1])
  // Outputs:
  //   h  hue value (degrees: [0,360])
  //   s  saturation value ([0,1])
  //   v  value value ([0,1])
  template <typename R,typename H>
  void rgb_to_hsv(const R * rgb, H * hsv);
};

#ifdef IGL_HEADER_ONLY
#  include "rgb_to_hsv.cpp"
#endif

#endif

