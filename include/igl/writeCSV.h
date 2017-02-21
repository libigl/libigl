// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Gustavo Segovia <gustavo.segovia@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_WRITE_CSV_H
#define IGL_WRITE_CSV_H

#include "igl/igl_inline.h"
#include <Eigen/Core>
#include <string>

namespace igl 
{
  // Write a CSV file based on matrix values
  // Inputs:
  //   str  path to outputfile
  //   M  eigen double matrix
  // Returns true on success, false on error
  template <typename Scalar>
  IGL_INLINE bool writeCSV(
	const std::string str, 
    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& M);
}

#ifndef IGL_STATIC_LIBRARY
#  include "writeCSV.cpp"
#endif

#endif
