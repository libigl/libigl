// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "parula.h"


template <typename T>
IGL_INLINE void igl::parula(const T x, T * rgb)
{
  return parula(x,rgb[0],rgb[1],rgb[2]);
}

template <typename T>
IGL_INLINE void igl::parula(const T f, T & r, T & g, T & b)
{
  // clamp to (0,1)
  const float eff_f = (f>=1.?1.:(f<=0.?0.:f));
  // continuous index into array
  const float ff = eff_f*(PARULA_COLOR_MAP.rows()-1);
  size_t s = floor(ff);
  size_t d = ceil(ff);
  const float t = (s==d ? 0. : (ff-s)/float(d-s));

  assert(t>=0 && t<=1);
  const auto rgb_s = PARULA_COLOR_MAP.row(s);
  const auto rgb_d = PARULA_COLOR_MAP.row(d);
  const auto rgb = rgb_s + t*(rgb_d-rgb_s);
  r = rgb(0);
  g = rgb(1);
  b = rgb(2);
}


template <typename DerivedZ, typename DerivedC>
IGL_INLINE void igl::parula(
  const Eigen::PlainObjectBase<DerivedZ> & Z,
  const bool normalize,
  Eigen::PlainObjectBase<DerivedC> & C)
{
  const double min_z = (normalize?Z.minCoeff():0);
  const double max_z = (normalize?Z.maxCoeff():1);
  return parula(Z,min_z,max_z,C);
}
template <typename DerivedZ, typename DerivedC>
IGL_INLINE void igl::parula(
  const Eigen::PlainObjectBase<DerivedZ> & Z,
  const double min_z,
  const double max_z,
  Eigen::PlainObjectBase<DerivedC> & C)
{
  C.resize(Z.rows(),3);
  for(int r = 0;r<Z.rows();r++)
  {
    parula((-min_z+Z(r,0))/(max_z-min_z),C(r,0),C(r,1),C(r,2));
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::parula<double>(double, double*);
template void igl::parula<Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, bool, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
