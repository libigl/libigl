//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

#ifndef IGL_WRITEOFF_H
#define IGL_WRITEOFF_H
#include "igl_inline.h"

#include <Eigen/Core>
#include <string>

namespace igl 
{
    IGL_INLINE bool writeOFF(const std::string fname, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
}

#ifdef IGL_HEADER_ONLY
#  include "writeOFF.cpp"
#endif

#endif
