//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

// History:
//  return type changed from void to bool  Alec 20 Sept 2011

#ifndef IGL_WRITEOBJ_H
#define IGL_WRITEOBJ_H
#include "igl_inline.h"

#include <Eigen/Core>
#include <string>

namespace igl 
{
  // Write a mesh in an ascii obj file
  // Inputs:
  //   str  path to outputfile
  //   V  eigen double matrix #V by 3 (mesh vertices)
  //   F  eigen int matrix #F by 3 (mesh indices)
  // Returns true on success, false on error
  IGL_INLINE bool writeOBJ(const std::string str, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
}

#ifdef IGL_HEADER_ONLY
#  include "writeOBJ.cpp"
#endif

#endif
