// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
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
