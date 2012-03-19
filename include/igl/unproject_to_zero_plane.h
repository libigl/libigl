#ifndef IGL_UNPROJECT_TO_ZERO_PLANE_H
#define IGL_UNPROJECT_TO_ZERO_PLANE_H
#include "igl_inline.h"
namespace igl
{
  // Wrapper for gluUnproject that uses the current GL_MODELVIEW_MATRIX,
  // GL_PROJECTION_MATRIX, and GL_VIEWPORT to unproject a screen postion
  // (winX,winY) to a 3d location at same depth as the current origin.
  // Inputs:
  //   win*  screen space x, y, and z coordinates respectively
  // Outputs:
  //   obj*  pointers to 3D objects' x, y, and z coordinates respectively
  // Returns return value of gluUnProject call
  IGL_INLINE int unproject_to_zero_plane(
    const double winX,
    const double winY,
    double* objX,
    double* objY,
    double* objZ);
}

#ifdef IGL_HEADER_ONLY
#  include "unproject_to_zero_plane.cpp"
#endif

#endif

