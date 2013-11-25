// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_CONCAT_H
#define IGL_CONCAT_H
#include "igl_inline.h"
#include <Eigen/Dense>
namespace igl
{
  // Concatenates dense matrices
  // Templates:
  //   T  should be a eigen matrix primitive type like int or double
  // Inputs:
  //   A first matrix
  //   B second matrix
  //   horiz if true, matrices are concatenated horizontally
  // Output:
  //   O if horiz = false return [A;B] else [A,B]
  template <typename T>
  IGL_INLINE void concat(
                     const T A, 
                     const T B,
                     const bool horiz,                 
                     T& O);
  
  template <typename T>
  IGL_INLINE T concat(
                  const T A, 
                  const T B,
                  bool horiz = false
                  );

}

#ifdef IGL_HEADER_ONLY
#  include "concat.cpp"
#endif

#endif
