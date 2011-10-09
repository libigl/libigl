#ifndef IGL_PROJECT_H
#define IGL_PROJECT_H
namespace igl
{
  // Wrapper for gluProject that uses the current GL_MODELVIEW_MATRIX,
  // GL_PROJECTION_MATRIX, and GL_VIEWPORT
  // Inputs:
  //   obj*  3D objects' x, y, and z coordinates respectively
  // Outputs:
  //   win*  pointers to screen space x, y, and z coordinates respectively
  // Returns return value of gluProject call
  inline int project(
    const double objX,
    const double objY,
    const double objZ,
    double* winX,
    double* winY,
    double* winZ);
}

// Implementation

#ifdef __APPLE__
# include <OpenGL/gl.h>
# include <OpenGL/glu.h>
#else
# include <GL/gl.h>
# include <GL/glu.h>
#endif

inline int igl::project(
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
#endif
