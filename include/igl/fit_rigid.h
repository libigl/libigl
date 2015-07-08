// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_FIT_RIGID_H
#define IGL_FIT_RIGID_H

#if defined(_WIN32)
  #pragma message("Deprecated. Use igl/procrustes.h instead")
#else
  #warning "Deprecated. Use igl/procrustes.h instead"
#endif

#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace igl
{
  // Deprecated: please use procrustes(...) instead.
  // Fit a rigid 
  IGL_INLINE void fit_rigid(
    const Eigen::MatrixXd & A,
    const Eigen::MatrixXd & B,
    Eigen::Rotation2Dd & R,
    Eigen::RowVector2d & t);
}

#ifndef IGL_STATIC_LIBRARY
#  include "fit_rigid.cpp"
#endif

#endif

