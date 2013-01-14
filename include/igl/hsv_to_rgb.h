#ifndef IGL_HSV_TO_RGB_H
#define IGL_HSV_TO_RGB_H
#include "igl_inline.h"
namespace igl
{
  // Convert RGB to HSV
  //
  // Inputs:
  //   h  hue value (degrees: [0,360])
  //   s  saturation value ([0,1])
  //   v  value value ([0,1])
  // Outputs:
  //   r  red value ([0,1]) 
  //   g  green value ([0,1])
  //   b  blue value ([0,1])
  template <typename T>
  void hsv_to_rgb(const T * hsv, T * rgb);
  template <typename T>
  void hsv_to_rgb( 
    const T & h, const T & s, const T & v, 
    T & r, T & g, T & b);
};

#ifdef IGL_HEADER_ONLY
#  include "hsv_to_rgb.cpp"
#endif

#endif

