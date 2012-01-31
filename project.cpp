#include "project.h"

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#else
#  ifdef _WIN32
#    define NOMINMAX
#    include <Windows.h>
#    undef NOMINMAX
#  endif
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif

IGL_INLINE int igl::project(
  const double objX,
  const double objY,
  const double objZ,
  double* winX,
  double* winY,
  double* winZ)
{
  // Put model, projection, and viewport matrices into double arrays
  double MV[16];
  double P[16];
  int VP[4];
  glGetDoublev(GL_MODELVIEW_MATRIX,  MV);
  glGetDoublev(GL_PROJECTION_MATRIX, P);
  glGetIntegerv(GL_VIEWPORT, VP);
  return gluProject(objX,objY,objZ,MV,P,VP,winX,winY,winZ);
}
