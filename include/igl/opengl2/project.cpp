// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "project.h"
#include "../project.h"
#include "../opengl/report_gl_error.h"
#include "../opengl/OpenGL_convenience.h"
#include <iostream>

IGL_INLINE int igl::opengl2::project(
  const double objX,
  const double objY,
  const double objZ,
  double* winX,
  double* winY,
  double* winZ)
{
  using namespace std;
#ifdef EXTREME_VERBOSE
  cout<<"project();"<<endl;
#endif
  // Put model, projection, and viewport matrices into double arrays
  double MV[16];
  double P[16];
  int VP[4];
  glGetDoublev(GL_MODELVIEW_MATRIX,  MV);

#ifdef EXTREME_VERBOSE
  cout<<"MV=["<<endl<<
    MV[0]<<" "<< MV[1]<<" "<< MV[2]<<" "<< MV[3]<<" "<<endl<<
    MV[4]<<" "<< MV[5]<<" "<< MV[6]<<" "<< MV[7]<<" "<<endl<<
    MV[8]<<" "<< MV[9]<<" "<< MV[10]<<" "<< MV[11]<<" "<<endl<<
    MV[12]<<" "<< MV[13]<<" "<< MV[14]<<" "<< MV[15]<<" "<<endl<<
    "];"<<endl;
#endif
#ifndef NDEBUG
  igl::opengl::report_gl_error();
#endif

  glGetDoublev(GL_PROJECTION_MATRIX, P);

#ifdef EXTREME_VERBOSE
  cout<<"P=["<<endl<<
    P[0]<<" "<< P[1]<<" "<< P[2]<<" "<< P[3]<<" "<<endl<<
    P[4]<<" "<< P[5]<<" "<< P[6]<<" "<< P[7]<<" "<<endl<<
    P[8]<<" "<< P[9]<<" "<< P[10]<<" "<< P[11]<<" "<<endl<<
    P[12]<<" "<< P[13]<<" "<< P[14]<<" "<< P[15]<<" "<<endl<<
    "];"<<endl;
#endif
#ifndef NDEBUG
  igl::opengl::report_gl_error();
#endif

  glGetIntegerv(GL_VIEWPORT, VP);

#ifdef EXTREME_VERBOSE
  cout<<"VP=["<<endl<<
    VP[0]<<" "<< VP[1]<<" "<< VP[2]<<" "<< VP[3]<<" "<<endl<<
    "];"<<endl;
#endif
#ifndef NDEBUG
  igl::opengl::report_gl_error();
#endif

#ifdef EXTREME_VERBOSE
  cout<<"obj=["<<endl<<
    objX<<" "<< objY<<" "<< objZ<<endl<<
    "];"<<endl;
#endif


  int ret = gluProject(objX,objY,objZ,MV,P,VP,winX,winY,winZ);

#ifdef EXTREME_VERBOSE
  cout<<"win=["<<endl<<
    *winX<<" "<< *winY<<" "<< *winZ<<endl<<
    "];"<<endl;
#endif
  return ret;
}

template <typename Derivedobj, typename Derivedwin>
IGL_INLINE int igl::opengl2::project(
  const Eigen::PlainObjectBase<Derivedobj> & obj,
  Eigen::PlainObjectBase<Derivedwin> & win)
{
  assert(obj.size() >= 3);
  Eigen::Vector3d dobj(obj(0),obj(1),obj(2));
  Eigen::Vector3d dwin;
  int ret = igl::opengl2::project(dobj(0),dobj(1),dobj(2),
      &dwin.data()[0],
      &dwin.data()[1],
      &dwin.data()[2]);
  win(0) = dwin(0);
  win(1) = dwin(1);
  win(2) = dwin(2);
  return ret;
}

template <typename Derivedobj>
IGL_INLINE Eigen::PlainObjectBase<Derivedobj> igl::opengl2::project(
  const Eigen::PlainObjectBase<Derivedobj> & obj)
{
  Eigen::PlainObjectBase<Derivedobj> win;
  igl::opengl2::project(obj,win);
  return win;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template int igl::opengl2::project<Eigen::Matrix<double, 3, 1, 0, 3, 1>, Eigen::Matrix<double, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> >&);
template Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > igl::opengl2::project<Eigen::Matrix<double, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&);
template int igl::opengl2::project<Eigen::Matrix<float, 3, 1, 0, 3, 1>, Eigen::Matrix<float, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> >&);
template Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> > igl::opengl2::project<Eigen::Matrix<float, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<float, 3, 1, 0, 3, 1> > const&);
template Eigen::PlainObjectBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > igl::opengl2::project<Eigen::Matrix<double, 1, -1, 1, 1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&);
template int igl::opengl2::project<Eigen::Matrix<double, 1, 3, 1, 1, 3>, Eigen::Matrix<double, 1, 3, 1, 1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&);
template Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> > igl::opengl2::project<Eigen::Matrix<double, 1, 3, 1, 1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> > const&);
template Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > igl::opengl2::project<Eigen::Matrix<double, 1, 2, 1, 1, 2> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> > const&);
template Eigen::PlainObjectBase<Eigen::Matrix<double, 2, 1, 0, 2, 1> > igl::opengl2::project<Eigen::Matrix<double, 2, 1, 0, 2, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 2, 1, 0, 2, 1> > const&);
#endif

