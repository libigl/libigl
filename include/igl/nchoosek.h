// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_NCHOOSEK
#define IGL_NCHOOSEK
#include "igl_inline.h"
#include <vector>

#include <Eigen/Core>

namespace igl {
  IGL_INLINE void nchoosek(int offset,
                           int k,
                           int N,
                           std::vector<std::vector<int> > &allCombs);
}


#ifndef IGL_STATIC_LIBRARY
#include "nchoosek.cpp"
#endif


#endif /* defined(IGL_NCHOOSEK) */
