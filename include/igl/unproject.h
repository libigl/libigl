// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_UNPROJECT_H
#define IGL_UNPROJECT_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Eigen reimplementation of gluUnproject
  // Inputs:
  //   win  screen space x, y, and z coordinates
  // Returns:
  //   the unprojected x, y, and z coordinates
  // Returns return value of gluUnProject call
  template <typename Scalar>
  IGL_INLINE Eigen::Matrix<Scalar,3,1> unproject(
    const    Eigen::Matrix<Scalar,3,1>&  win,
    const    Eigen::Matrix<Scalar,4,4>& model,
    const    Eigen::Matrix<Scalar,4,4>& proj,
    const    Eigen::Matrix<Scalar,4,1>&  viewport);
}


#ifndef IGL_STATIC_LIBRARY
#  include "unproject.cpp"
#endif

#endif
