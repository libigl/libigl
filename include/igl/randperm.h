// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_RANDPERM_H
#define IGL_RANDPERM_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <functional>
namespace igl
{
  // Like matlab's randperm(n) but minus 1
  //
  // Inputs:
  //   n  number of elements
  //   rng_min  the minimum value of rng()
  //   rng_max  the maximum value of rng()
  //   rng random number generator. When not given the value,
  //       randperm will use default random number generator std::rand()
  // Outputs:
  //   I  n list of rand permutation of 0:n-1
  template <typename DerivedI>
  IGL_INLINE void randperm(
    const int n,
    Eigen::PlainObjectBase<DerivedI> & I,
    const int64_t rng_min=0,
    const int64_t rng_max=RAND_MAX,
    const std::function<int64_t()> &rng=nullptr);
}
#ifndef IGL_STATIC_LIBRARY
#  include "randperm.cpp"
#endif
#endif
