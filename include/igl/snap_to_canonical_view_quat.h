// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SNAP_TO_CANONICAL_VIEW_QUAT_H
#define IGL_SNAP_TO_CANONICAL_VIEW_QUAT_H
#include "igl_inline.h"
#include <Eigen/Geometry>
// A Quaternion, q, is defined here as an arrays of four scalars (x,y,z,w),
// such that q = x*i + y*j + z*k + w
namespace igl
{
  // Snap a given quaternion to the "canonical quaternion" rotations.
  // Inputs:
  //   q  input quaternion
  //   threshold  threshold between 0 and 1, where 0 means
  template <typename Q_type>
  IGL_INLINE bool snap_to_canonical_view_quat(
    const Q_type q[4],
    const Q_type threshold,
    Q_type s[4]);
  IGL_INLINE bool snap_to_canonical_view_quat(
    const Eigen::Quaterniond & q,
    const double threshold,
    Eigen::Quaterniond & s);
}

#ifdef IGL_HEADER_ONLY
#  include "snap_to_canonical_view_quat.cpp"
#endif

#endif
