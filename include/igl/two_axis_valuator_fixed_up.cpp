// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "two_axis_valuator_fixed_up.h"
#include "PI.h"

IGL_INLINE void igl::two_axis_valuator_fixed_up(
  const int w,
  const int h,
  const double speed,
  const Eigen::Quaterniond & down_quat,
  const int down_x,
  const int down_y,
  const int mouse_x,
  const int mouse_y,
  Eigen::Quaterniond & quat)
{
  using namespace Eigen;
  Vector3d axis(0,1,0);
  quat = down_quat * 
    Quaterniond(
      AngleAxisd(
        PI*((double)(mouse_x-down_x))/(double)w*speed/2.0,
        axis.normalized()));
  quat.normalize();
  {
    Vector3d axis(1,0,0);
    if(axis.norm() != 0)
    {
      quat = 
        Quaterniond(
          AngleAxisd(
            PI*(mouse_y-down_y)/(double)h*speed/2.0,
            axis.normalized())) * quat;
      quat.normalize();
    }
  }
}

