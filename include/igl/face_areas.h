// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_FACE_AREAS_H
#define IGL_FACE_AREAS_H
#include "igl_inline.h"

#include <Eigen/Dense>

namespace igl
{
  // Constructs a list of face areas of faces opposite each index in a tet list
  //
  // Templates:
  //   DerivedV derived from vertex positions matrix type: i.e. MatrixXd
  //   DerivedT derived from face indices matrix type: i.e. MatrixXi
  //   DerivedL derived from edge lengths matrix type: i.e. MatrixXd
  // Inputs:
  //   V  eigen matrix #V by 3
  //   T  #T by 3 list of mesh faces (must be triangles)
  // Outputs:
  //   E #E by 2 list of edges in no particular order
  //
  // See also: adjacency_matrix
  template <typename DerivedV, typename DerivedT, typename DerivedA>
  IGL_INLINE void face_areas(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedT>& T,
    Eigen::PlainObjectBase<DerivedA>& A);
  template <typename DerivedL, typename DerivedA>
  IGL_INLINE void face_areas(
    const Eigen::PlainObjectBase<DerivedL>& L,
    Eigen::PlainObjectBase<DerivedA>& A);
}

#ifndef IGL_STATIC_LIBRARY
#  include "face_areas.cpp"
#endif

#endif


