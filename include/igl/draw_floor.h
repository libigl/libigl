#ifndef IGL_DRAW_FLOOR_H
#define IGL_DRAW_FLOOR_H
#ifndef IGL_NO_OPENGL
#include "igl_inline.h"
namespace igl
{
  // Draw a checkerboard floor aligned with current (X,Z) plane using OpenGL
  // calls.
  //
  // Use glPushMatrix(), glScaled(), glTranslated() to arrange the floor.
  // 
  // Inputs:
  //   colorA  float4 color
  //   colorB  float4 color
  IGL_INLINE void draw_floor(const float * colorA, const float * colorB);
  // Wrapper with default colors
  IGL_INLINE void draw_floor();
}
#ifdef IGL_HEADER_ONLY
#  include "draw_floor.cpp"
#endif
#endif
#endif
