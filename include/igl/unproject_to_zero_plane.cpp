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

template <typename Derivedwin, typename Derivedobj>
IGL_INLINE int igl::unproject_to_zero_plane(
  const Eigen::PlainObjectBase<Derivedwin> & win,
  Eigen::PlainObjectBase<Derivedobj> & obj)
{
  return unproject_to_zero_plane(win(0),win(1),
      &obj.data()[0],
      &obj.data()[1],
      &obj.data()[2]);
}

#ifndef IGL_HEADER_ONLY
// Explicit template instanciation
template int igl::unproject_to_zero_plane<Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::Matrix<double, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> >&);
#endif

#endif
