//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

#ifndef IGL_WRITE_H
#define IGL_WRITE_H
#include "igl_inline.h"

#include <Eigen/Core>
#include <string>

namespace igl 
{
  // write mesh to an ascii file with automatic detection of file format. supported: obj, off)
  // Known Bugs:
  //  Does not correctly find file extensions: myfile.foo.off 
  IGL_INLINE bool write(
    const std::string str, 
    const Eigen::MatrixXd& V, 
    const Eigen::MatrixXi& F);
}

#ifdef IGL_HEADER_ONLY
#  include "write.cpp"
#endif

#endif
