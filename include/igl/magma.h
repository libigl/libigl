// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2017 Joe Graus <jgraus@gmu.edu>, Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_MAGMA_H
#define IGL_MAGMA_H
#include "igl_inline.h"

#include <Eigen/Dense>

namespace igl {
  // Magma colormap -- a perceptually uniform sequential colormap
  //
  // Inputs:
  //   m  number of colors 
  // Outputs:
  //   J  m by list of RGB colors between 0 and 1
  //
  // Wrapper for directly computing [r,g,b] values for a given factor f between
  // 0 and 1
  //
  // Inputs:
  //   f  factor determining color value as if 0 was min and 1 was max
  // Outputs:
  //   r  red value
  //   g  green value
  //   b  blue value
  template <typename T>
  IGL_INLINE void magma(const T f, T * rgb);
  template <typename T>
  IGL_INLINE void magma(const T f, T & r, T & g, T & b);
  // Inputs:
  //   Z  #Z list of factors
  //   normalize  whether to normalize Z to be tightly between [0,1]
  // Outputs:
  //   C  #C by 3 list of rgb colors
  template <typename DerivedZ, typename DerivedC>
  IGL_INLINE void magma(
    const Eigen::PlainObjectBase<DerivedZ> & Z,
    const bool normalize,
    Eigen::PlainObjectBase<DerivedC> & C);
  // Inputs:
  //   min_z  value at black
  //   max_z  value at yellow
  template <typename DerivedZ, typename DerivedC>
  IGL_INLINE void magma(
    const Eigen::PlainObjectBase<DerivedZ> & Z,
    const double min_Z,
    const double max_Z,
    Eigen::PlainObjectBase<DerivedC> & C);
};

#ifndef IGL_STATIC_LIBRARY
#  include "magma.cpp"
#endif

#endif
