// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
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

template <typename Derivedwin, typename Derivedobj>
IGL_INLINE int igl::unproject(
  const Eigen::PlainObjectBase<Derivedwin> & win,
  Eigen::PlainObjectBase<Derivedobj> & obj)
{
  Eigen::Vector3d dwin(win(0),win(1),win(2));
  Eigen::Vector3d dobj;
  int ret = unproject(dwin(0),dwin(1),dwin(2),
      &dobj.data()[0],
      &dobj.data()[1],
      &dobj.data()[2]);
  obj(0) = dobj(0);
  obj(1) = dobj(1);
  obj(2) = dobj(2);
  return ret;
}

template <typename Derivedwin>
IGL_INLINE Eigen::PlainObjectBase<Derivedwin> igl::unproject(
  const Eigen::PlainObjectBase<Derivedwin> & win)
{
  Eigen::PlainObjectBase<Derivedwin> obj;
  unproject(win,obj);
  return obj;
}


#ifndef IGL_HEADER_ONLY
// Explicit template instanciation
template int igl::unproject<Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::Matrix<double, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> >&);
template int igl::unproject<Eigen::Matrix<float, 3, 1, 0, 3, 1>, Eigen::Matrix<float, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> >&);
template Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> > igl::unproject<Eigen::Matrix<float, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> > const&);
template int igl::unproject<Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::Matrix<float, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> >&);
template Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > igl::unproject<Eigen::Matrix<double, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&);
template int igl::unproject<Eigen::Matrix<double, 1, 3, 1, 1, 3>, Eigen::Matrix<double, 1, 3, 1, 1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&);
#endif

#endif
