// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2017 Joe Graus <jgraus@gmu.edu>, Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_COLORMAP_H
#define IGL_COLORMAP_H
#include "igl_inline.h"

#include <Eigen/Dense>

namespace igl {
  // Colormap selector -- an interface to supported colormaps within igl
  //
  // Inputs:
  //   m  number of colors 
  //   cm colormap enum
  // Outputs:
  //   J  m by list of RGB colors between 0 and 1
  //
  // Wrapper for directly computing [r,g,b] values of the selected colormap for a given factor f between
  // 0 and 1
  //
  // Inputs:
  //   f  factor determining color value as if 0 was min and 1 was max
  //   c  colormap enum
  // Outputs:
  //   r  red value
  //   g  green value
  //   b  blue value

  enum ColorMapType
  {
	CM_INFERNO = 0,
	CM_JET = 1,
	CM_MAGMA = 2,
	CM_PARULA = 3,
	CM_PLASMA = 4,
	CM_VIRIDIS = 5
  };

  template <typename T>
  IGL_INLINE void colormap(const T f, T * rgb, ColorMapType cm);
  template <typename T>
  IGL_INLINE void colormap(const T f, T & r, T & g, T & b, ColorMapType cm);
  // Inputs:
  //   Z  #Z list of factors
  //   normalize  whether to normalize Z to be tightly between [0,1]
  //   cm selected colormap palette to interpolate from
  // Outputs:
  //   C  #C by 3 list of rgb colors
  template <typename DerivedZ, typename DerivedC>
  IGL_INLINE void colormap(
    const Eigen::PlainObjectBase<DerivedZ> & Z,
    const bool normalize,
    Eigen::PlainObjectBase<DerivedC> & C,
	ColorMapType cm);
  // Inputs:
  //   min_z  value at black
  //   max_z  value at yellow
  template <typename DerivedZ, typename DerivedC>
  IGL_INLINE void colormap(
    const Eigen::PlainObjectBase<DerivedZ> & Z,
    const double min_Z,
    const double max_Z,
    Eigen::PlainObjectBase<DerivedC> & C,
	ColorMapType cm);
};

#ifndef IGL_STATIC_LIBRARY
#  include "colormap.cpp"
#endif

#endif
