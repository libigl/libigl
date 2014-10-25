// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "draw_point.h"
#ifndef IGL_NO_OPENGL

// Implementation
#include "OpenGL_convenience.h"

#include <cassert>
#include <cmath>

IGL_INLINE void igl::draw_point(
  const double x,
  const double y,
  const double z,
  const double requested_r,
  const bool selected)
{
  // Push GL settings
  //GLboolean old_depth_test;
  //glGetBooleanv(GL_DEPTH_TEST,&old_depth_test);
  GLboolean old_lighting;
  glGetBooleanv(GL_LIGHTING,&old_lighting);
  glEnable( GL_POINT_SMOOTH );

  float f;
  glGetFloatv(GL_POINT_SIZE_MAX,&f);
  // THIS IS OVERZEALOUS on Mac OS X: OpenGL reports a smaller point size than
  // possible.
  //assert(requested_r<=0.5*f);
  double r = (requested_r<0.5*f?requested_r:0.5*f);

  //glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);

  // get current color
  float color[4];
  glGetFloatv(GL_CURRENT_COLOR,color);

  double outline_size = (r>7 ? sqrt(r/7.0) : 1.0);

  // White outline
  glColor4f(1,1,1,color[3]);
  glPointSize(2*r);
  glBegin(GL_POINTS);
  glVertex3d(x,y,z);
  glEnd();
  // Black outline
  glColor4f(0,0,0,color[3]);
  glPointSize(2*r-2*outline_size);
  glBegin(GL_POINTS);
  glVertex3d(x,y,z);
  glEnd();
  
  // Foreground
  glColor4fv(color);
  glPointSize(2*r-4*outline_size);
  glBegin(GL_POINTS);
  glVertex3d(x,y,z);
  glEnd();
  // Selection inner circle
  if(selected)
  {
    glColor4f(0,0,0,color[3]);
    double selected_size = 2*r-7*outline_size;
    selected_size = (selected_size>3?selected_size:3);
    glPointSize(selected_size);
    glBegin(GL_POINTS);
    glVertex3d(x,y,z);
    glEnd();
  }

  // reset color
  glColor4fv(color);

  // Pop GL settings
  if(old_lighting) glEnable(GL_LIGHTING);
  //if(old_depth_test) glEnable(GL_DEPTH_TEST);
}

template <typename DerivedP>
IGL_INLINE void igl::draw_point(
  const Eigen::PlainObjectBase<DerivedP> & P,
  const double requested_r,
  const bool selected)
{
  switch(P.size())
  {
    case 2:
      return draw_point(P(0),P(1),0,requested_r,selected);
    default:
      return draw_point(P(0),P(1),P(2),requested_r,selected);
  }
}

#ifdef IGL_STATIC_LIBRARY
template void igl::draw_point<Eigen::Matrix<double, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&, double, bool);
template void igl::draw_point<Eigen::Matrix<double, 2, 1, 0, 2, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 2, 1, 0, 2, 1> > const&, double, bool); 
#endif

#endif
