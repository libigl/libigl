// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2024 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/
#include "find_intersections.h"
#include "find_self_intersections.h"

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedIF>
IGL_INLINE bool igl::predicates::find_self_intersections(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  const bool first_only,
  Eigen::PlainObjectBase<DerivedIF> & IF)
{
  // This is really just a wrapper around fast_find_intersections which will
  // internally detect that V,F are the second set
  return find_intersections( V,F,V,F,first_only,IF);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template bool igl::predicates::find_self_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, bool, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
