//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

// History:
//  return type changed from void to bool  Alec 18 Sept 2011

#ifndef IGL_READOFF_H
#define IGL_READOFF_H
#include "igl_inline.h"

#include <Eigen/Core>
#include <string>

namespace igl 
{
  // read mesh from a ascii off file
  // Inputs:
  //   str  path to .off file
  // Outputs:
  //   V  eigen double matrix #V by 3
  //   F  eigen int matrix #F by 3
  IGL_INLINE bool readOFF (const std::string meshfile, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
}

#ifdef IGL_HEADER_ONLY
#  include "readOFF.cpp"
#endif

#endif
