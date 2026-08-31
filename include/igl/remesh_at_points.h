// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_REMESH_AT_POINTS_H
#define IGL_REMESH_AT_POINTS_H
#include "igl_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace igl 
{
  ///   @param[in] V  #V by dim list of mesh vertex positions
  ///   @param[in] F  #F by 3 list of triangle indices into rows of V
  ///   @param[in] B  #B by 3 list of barycentric coordinates, ith row are
  ///   coordinates of ith sampled point in face FI(i)
  ///   @param[in] FI  #B list of indices into F
  ///   @param[out] VV  #VV by dim list of mesh vertex positions (top #V rows is
  ///   always V; bottom #B rows are new vertices)
  ///   @param[out] FF  #FF by 3 list of triangle indices into rows of VV
  ///   @param[out] J  #FF list of indices into F 
  ///   @param[out] K  #B list of indices into VV
  template <
    typename DerivedV, 
    typename DerivedF, 
    typename DerivedB,
    typename DerivedFI,
    typename DerivedVV,
    typename DerivedFF,
    typename DerivedJ,
    typename DerivedK
    >
  IGL_INLINE void remesh_at_points(
    const Eigen::MatrixBase<DerivedV> & V, 
    const Eigen::MatrixBase<DerivedF> & F, 
    const Eigen::MatrixBase<DerivedB> & B,
    const Eigen::MatrixBase<DerivedFI> & FI,
    Eigen::PlainObjectBase<DerivedVV> & VV,
    Eigen::PlainObjectBase<DerivedFF> & FF,
    Eigen::PlainObjectBase<DerivedJ> & J,
    Eigen::PlainObjectBase<DerivedK> & K);
}

#ifndef IGL_STATIC_LIBRARY
#  include "remesh_at_points.cpp"
#endif

#endif

