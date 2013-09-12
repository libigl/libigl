#include "unproject.h"
#ifndef IGL_NO_OPENGL

#include "OpenGL_convenience.h"

IGL_INLINE int igl::unproject(
  const double winX,
  const double winY,
  const double winZ,
  double* objX,
  double* objY,
  double* objZ)
{
  // Put model, projection, and viewport matrices into double arrays
  double MV[16];
  double P[16];
  int VP[4];
  glGetDoublev(GL_MODELVIEW_MATRIX,  MV);
  glGetDoublev(GL_PROJECTION_MATRIX, P);
  glGetIntegerv(GL_VIEWPORT, VP);
  return gluUnProject(winX,winY,winZ,MV,P,VP,objX,objY,objZ);
}
#endif
