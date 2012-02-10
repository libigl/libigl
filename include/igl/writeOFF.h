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
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE bool writeOFF(
                        const std::string str,
                        const Eigen::PlainObjectBase<DerivedV>& V,
                        const Eigen::PlainObjectBase<DerivedF>& F);
}

#ifdef IGL_HEADER_ONLY
#  include "writeOFF.cpp"
#endif

#endif
