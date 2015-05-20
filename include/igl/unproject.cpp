// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "unproject.h"
#ifndef IGL_NO_OPENGL
#ifndef IGL_OPENGL_4

#include <Eigen/Dense>
#include <Eigen/LU>
#include "OpenGL_convenience.h"

IGL_INLINE void igl::unproject(
  const double winX,
  const double winY,
  const double winZ,
  double* objX,
  double* objY,
  double* objZ)
{
  Eigen::Vector3d obj;
  igl::unproject(Eigen::Vector3d(winX,winY,winZ),obj);
  *objX = obj(0);
  *objY = obj(1);
  *objZ = obj(2);
}

template <typename Derivedwin, typename Derivedobj>
IGL_INLINE void igl::unproject(
  const Eigen::PlainObjectBase<Derivedwin> & win,
  Eigen::PlainObjectBase<Derivedobj> & obj)
{
  obj = igl::unproject(win).template cast<typename Derivedobj::Scalar>();
}

template <typename Derivedwin>
IGL_INLINE Eigen::PlainObjectBase<Derivedwin> igl::unproject(
  const Eigen::PlainObjectBase<Derivedwin> & win)
{
  using namespace Eigen;
  typedef typename Derivedwin::Scalar Scalar;
  Matrix4d MV,P;
  Vector4i VPi;
  Vector4d VPd;
  glGetDoublev(GL_MODELVIEW_MATRIX,MV.data());
  glGetDoublev(GL_PROJECTION_MATRIX,P.data());
  glGetIntegerv(GL_VIEWPORT,VPi.data());
  VPd = VPi.cast<double>();
  Vector3d wind = win.template cast<double>();
  Vector3d objd = igl::unproject(wind,MV,P,VPd);
  return objd.template cast<Scalar>();
}

#endif
#endif


template <typename Scalar>
IGL_INLINE Eigen::Matrix<Scalar,3,1> igl::unproject(
  const    Eigen::Matrix<Scalar,3,1>&  win,
  const    Eigen::Matrix<Scalar,4,4>& model,
  const    Eigen::Matrix<Scalar,4,4>& proj,
  const    Eigen::Matrix<Scalar,4,1>&  viewport)
{
  Eigen::Matrix<Scalar,4,4> Inverse = (proj * model).inverse();

  Eigen::Matrix<Scalar,4,1> tmp;
  tmp << win, 1;
  tmp(0) = (tmp(0) - viewport(0)) / viewport(2);
  tmp(1) = (tmp(1) - viewport(1)) / viewport(3);
  tmp = tmp.array() * 2.0f - 1.0f;

  Eigen::Matrix<Scalar,4,1> obj = Inverse * tmp;
  obj /= obj(3);

  return obj.head(3);
}

#ifdef IGL_STATIC_LIBRARY

#ifndef IGL_NO_OPENGL
#ifndef IGL_OPENGL_4
// Explicit template instanciation
template void igl::unproject<Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::Matrix<double, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> >&);
template void igl::unproject<Eigen::Matrix<float, 3, 1, 0, 3, 1>, Eigen::Matrix<float, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> >&);
template Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> > igl::unproject<Eigen::Matrix<float, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> > const&);
template void igl::unproject<Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::Matrix<float, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> >&);
template Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > igl::unproject<Eigen::Matrix<double, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&);
template void igl::unproject<Eigen::Matrix<double, 1, 3, 1, 1, 3>, Eigen::Matrix<double, 1, 3, 1, 1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&);
template Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> > igl::unproject<Eigen::Matrix<double, 1, 3, 1, 1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> > const&);
template Eigen::Matrix<float, 3, 1, 0, 3, 1> igl::unproject<float>(Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, Eigen::Matrix<float, 4, 1, 0, 4, 1> const&);
#endif
template Eigen::Matrix<double, 3, 1, 0, 3, 1> igl::unproject<double>(Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 4, 4, 0, 4, 4> const&, Eigen::Matrix<double, 4, 4, 0, 4, 4> const&, Eigen::Matrix<double, 4, 1, 0, 4, 1> const&);
#endif

#endif
