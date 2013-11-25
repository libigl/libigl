// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PER_VERTEX_NORMALS_H
#define IGL_PER_VERTEX_NORMALS_H
#include "igl_inline.h"
#include <Eigen/Core>
// Note: So for this only computes normals per vertex as uniformly weighted
// averages of incident triangle normals. It would be nice to support more or
// all of the methods here:
// "A comparison of algorithms for vertex normal computation"
namespace igl
{
  // Compute vertex normals via vertex position list, face list
  // Inputs:
  //   V  #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F  #F by 3 eigne Matrix of face (triangle) indices
  // Output:
  //   N  #V by 3 eigen Matrix of mesh vertex 3D normals
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE void per_vertex_normals(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedV> & N);

  template <typename DerivedV, typename DerivedF>
  IGL_INLINE void per_vertex_normals(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F,
    const Eigen::PlainObjectBase<DerivedV>& FN,
    Eigen::PlainObjectBase<DerivedV> & N);

}

#ifdef IGL_HEADER_ONLY
#  include "per_vertex_normals.cpp"
#endif

#endif
