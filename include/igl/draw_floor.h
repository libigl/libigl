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
  //
  // Example:
  //   // Draw a nice floor
  //   glPushMatrix();
  //   glCullFace(GL_BACK);
  //   glEnable(GL_CULL_FACE);
  //   glEnable(GL_LIGHTING);
  //   glTranslated(0,-1,0);
  //   if(project(Vector3d(0,0,0))(2) - project(Vector3d(0,1,0))(2) > -FLOAT_EPS)
  //   {
  //     draw_floor_outline();
  //   }
  //   draw_floor();
  //   glPopMatrix();
  //   glDisable(GL_CULL_FACE);
  //
  IGL_INLINE void draw_floor(const float * colorA, const float * colorB);
  // Wrapper with default colors
  IGL_INLINE void draw_floor();
  IGL_INLINE void draw_floor_outline(const float * colorA, const float * colorB);
  // Wrapper with default colors
  IGL_INLINE void draw_floor_outline();
}
#ifdef IGL_HEADER_ONLY
#  include "draw_floor.cpp"
#endif
#endif
#endif
