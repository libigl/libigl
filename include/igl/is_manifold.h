//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

#ifndef IGL_IS_MANIFOLD_H
#define IGL_IS_MANIFOLD_H
#include "igl_inline.h"

#include <Eigen/Core>
#include <vector>

namespace igl 
{
  // check if the mesh is edge-manifold
  //
  // Not clear whether this returns true or false if the mesh is disc topology
  //
  // Known Bugs:
  //  Does not check for non-manifold vertices
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE bool is_manifold(const Eigen::PlainObjectBase<DerivedV>& V,
                              const Eigen::PlainObjectBase<DerivedF>& F);
}

#ifdef IGL_HEADER_ONLY
#  include "is_manifold.cpp"
#endif

#endif
