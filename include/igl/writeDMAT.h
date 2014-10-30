// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_WRITEDMAT_H
#define IGL_WRITEDMAT_H
#include "igl_inline.h"
// See writeDMAT.h for a description of the .dmat file type
#include <string>
#include <vector>
namespace igl
{
  // Write a matrix using ascii dmat file type
  //
  // Template:
  //   Mat  matrix type that supports .rows(), .cols(), operator(i,j)
  // Inputs:
  //   file_name  path to .dmat file
  //   W  eigen matrix containing to-be-written coefficients
  //   ascii  write ascii file {true}
  // Returns true on success, false on error
  //
  template <class Mat>
  IGL_INLINE bool writeDMAT(
    const std::string file_name, 
    const Mat & W,
    const bool ascii=true);
  template <typename Scalar>
  IGL_INLINE bool writeDMAT(
    const std::string file_name, 
    const std::vector<std::vector<Scalar> > W);
  template <typename Scalar>
  IGL_INLINE bool writeDMAT(
    const std::string file_name, 
    const std::vector<Scalar > W);
}

#ifndef IGL_STATIC_LIBRARY
#  include "writeDMAT.cpp"
#endif

#endif
