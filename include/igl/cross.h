// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_CROSS_H
#define IGL_CROSS_H
#include "igl_inline.h"
namespace igl
{
  // Computes out = cross(a,b)
  // Inputs:
  //   a  left 3d vector
  //   b  right 3d vector
  // Outputs:
  //   out  result 3d vector
  IGL_INLINE void cross(
    const double *a, 
    const double *b,
    double *out);
}

#ifdef IGL_HEADER_ONLY
#  include "cross.cpp"
#endif

#endif
