#include "draw_point.h"

// Implementation

#if __APPLE__
#  include <OpenGL/gl.h>
#elif defined(_WIN32)
#    define NOMINMAX
#    include <Windows.h>
#    undef NOMINMAX
#    include <GL/glew.h>
#    include <GL/gl.h>
#else
#  define GL_GLEXT_PROTOTYPES
#  include <GL/gl.h>
#  include <GL/glext.h>
#endif

#include <cassert>
#include <cmath>

IGL_INLINE void igl::draw_point(
  const double x,
  const double y,
  const double z,
  const double requested_r,
  const bool selected)
{
  // Push GL settings
  //GLboolean old_depth_test;
  //glGetBooleanv(GL_DEPTH_TEST,&old_depth_test);
  GLboolean old_lighting;
  glGetBooleanv(GL_LIGHTING,&old_lighting);

  float f;
  glGetFloatv(GL_POINT_SIZE_MAX,&f);
  assert(requested_r<=0.5*f);
  double r = (requested_r<0.5*f?requested_r:0.5*f);

  //glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);

  // get current color
  float color[4];
  glGetFloatv(GL_CURRENT_COLOR,color);

  double outline_size = (r>7 ? sqrt(r/7.0) : 1.0);

  // White outline
  glColor4f(1,1,1,color[3]);
  glPointSize(2*r);
  glBegin(GL_POINTS);
  glVertex3d(x,y,z);
  glEnd();
  // Black outline
  glColor4f(0,0,0,color[3]);
  glPointSize(2*r-2*outline_size);
  glBegin(GL_POINTS);
  glVertex3d(x,y,z);
  glEnd();
  
  // Foreground
  glColor4fv(color);
  glPointSize(2*r-4*outline_size);
  glBegin(GL_POINTS);
  glVertex3d(x,y,z);
  glEnd();
  // Selection inner circle
  if(selected)
  {
    glColor4f(0,0,0,color[3]);
    double selected_size = 2*r-7*outline_size;
    selected_size = (selected_size>3?selected_size:3);
    glPointSize(selected_size);
    glBegin(GL_POINTS);
    glVertex3d(x,y,z);
    glEnd();
  }

  // reset color
  glColor4fv(color);

  // Pop GL settings
  if(old_lighting) glEnable(GL_LIGHTING);
  //if(old_depth_test) glEnable(GL_DEPTH_TEST);
}

