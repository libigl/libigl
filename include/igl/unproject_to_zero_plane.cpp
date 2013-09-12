#include "unproject_to_zero_plane.h"
#ifndef IGL_NO_OPENGL

#include "OpenGL_convenience.h"

#include "project.h"
#include "unproject.h"

IGL_INLINE int igl::unproject_to_zero_plane(
  const double winX,
  const double winY,
  double* objX,
  double* objY,
  double* objZ)
{
  double winOrigin[3]; 
  igl::project(0,0,0,&winOrigin[0],&winOrigin[1],&winOrigin[2]);
  return igl::unproject(winX, winY, winOrigin[2], objX, objY, objZ);
}

#endif
