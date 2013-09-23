#ifndef IGL_UNPROJECT_H
#define IGL_UNPROJECT_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Wrapper for gluUnproject that uses the current GL_MODELVIEW_MATRIX,
  // GL_PROJECTION_MATRIX, and GL_VIEWPORT
  // Inputs:
  //   win*  screen space x, y, and z coordinates respectively
  // Outputs:
  //   obj*  pointers to 3D objects' x, y, and z coordinates respectively
  // Returns return value of gluUnProject call
  IGL_INLINE int unproject(
    const double winX,
    const double winY,
    const double winZ,
    double* objX,
    double* objY,
    double* objZ);
  template <typename Derivedwin, typename Derivedobj>
  IGL_INLINE int unproject(
    const Eigen::PlainObjectBase<Derivedwin> & win,
    Eigen::PlainObjectBase<Derivedobj> & obj);
}

#ifdef IGL_HEADER_ONLY
#  include "unproject.cpp"
#endif

#endif
