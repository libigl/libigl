#ifndef IGL_UNPROJECT_H
#define IGL_UNPROJECT_H
namespace igl
{
  // Wrapper for gluProject that uses the current GL_MODELVIEW_MATRIX,
  // GL_PROJECTION_MATRIX, and GL_VIEWPORT
  // Inputs:
  //   win*  screen space x, y, and z coordinates respectively
  // Outputs:
  //   obj*  pointers to 3D objects' x, y, and z coordinates respectively
  // Returns return value of gluUnProject call
  inline int unproject(
    const double winX,
    const double winY,
    const double winZ,
    double* objX,
    double* objY,
    double* objZ);
}

// Implementation

#ifdef __APPLE__
# include <OpenGL/gl.h>
# include <OpenGL/glu.h>
#else
#  ifdef _WIN32
#    define NOMINMAX
#    include <Windows.h>
#    undef NOMINMAX
#  endif
# include <GL/gl.h>
# include <GL/glu.h>
#endif

inline int igl::unproject(
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
