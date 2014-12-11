// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Stefan Brugger <stefanbrugger@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_OUTLINE_ORDERED_H
#define IGL_OUTLINE_ORDERED_H
#include "igl_inline.h"

#include <Eigen/Dense>
#include <vector>

namespace igl
{
  // Compute list of ordered boundary loops for a manifold mesh.
  //
  // Templates:
  //  Index  index type
  // Inputs:
  //   F  #V by dim list of mesh faces
  // Outputs:
  //   L  list of loops where L[i] = ordered list of boundary vertices in loop i
  //
  template <typename DerivedF, typename Index>
  IGL_INLINE void outline_ordered(
    const Eigen::PlainObjectBase<DerivedF> & F, 
    std::vector<std::vector<Index> >& L);
}

#ifndef IGL_STATIC_LIBRARY
#  include "outline_ordered.cpp"
#endif
#endif
