// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_ANGLES_H
#define IGL_ANGLES_H
#ifdef _WIN32
 #pragma message ( "Deprecated. Use igl/internal_angles.h instead" )
#else
 #warning "Deprecated. Use igl/internal_angles.h instead"
#endif


#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // ANGLES Compute angles for each corner of each triangle
  //
  // Inputs:
  //   V  #V by dim list of vertex positions
  //   F  #V by 3[4] list of triangle[quads] indices
  // Outputs:
  //   theta  #F by 3[4] list of angles for each corner (in radians)
  //
  template <
  typename DerivedV,
  typename DerivedF,
  typename Derivedtheta>
  IGL_INLINE void angles(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedF>& F,
  Eigen::PlainObjectBase<Derivedtheta>& theta);

}

#ifndef IGL_STATIC_LIBRARY
#  include "angles.cpp"
#endif

#endif
