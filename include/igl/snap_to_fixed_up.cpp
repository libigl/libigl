// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "snap_to_fixed_up.h"

IGL_INLINE void igl::snap_to_fixed_up(
  const Eigen::Quaterniond & q,
  Eigen::Quaterniond & s)
{
  using namespace Eigen;
  const Vector3d up = q.matrix() * Vector3d(0,1,0);
  Vector3d proj_up(0,up(1),up(2));
  if(proj_up.norm() == 0)
  {
    proj_up = Vector3d(0,1,0);
  }
  proj_up.normalize();
  Quaterniond dq;
  dq = Quaterniond::FromTwoVectors(up,proj_up);
  s = dq * q;
}
