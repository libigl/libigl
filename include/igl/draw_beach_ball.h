#ifndef IGL_DRAW_BEACH_BALL_H
#define IGL_DRAW_BEACH_BALL_H
#include "igl_inline.h"

namespace igl
{
  // Draw a beach ball icon/glyph (from AntTweakBar) at the current origin
  // according to the current orientation: ball has radius 0.75 and axis have
  // length 1.15
  IGL_INLINE void draw_beach_ball();
  #ifdef IGL_HEADER
  #  include "draw_beach_ball.cpp"
  #endif
}

#endif
