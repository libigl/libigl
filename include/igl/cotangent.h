// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_COTANGENT_H
#define IGL_COTANGENT_H
#include "igl_inline.h"
namespace igl
{
  // COTANGENT compute the cotangents of each angle in mesh (V,F)
  // 
  // Templates:
  //   MatV  vertex position matrix, e.g. Eigen::MatrixXd
  //   MatF  face index matrix, e.g. Eigen::MatrixXd
  //   MatC  cotangent weights matrix, e.g. Eigen::MatrixXd
  // Inputs:
  //   V  #V by dim list of rest domain positions
  //   F  #F by {3|4} list of {triangle|tetrahedra} indices into V
  // Outputs:
  //   C  #F by {3|6} list of cotangents corresponding angles
  //     for triangles, columns correspond to edges [1,2],[2,0],[0,1]
  //     for tets, columns correspond to edges [1,2],[2,0],[0,1],[3,0],[3,1],[3,2]
  //
  // Known bug:
  //   This computes 0.5*cotangent
  template <class MatV, class MatF, class MatC>
  IGL_INLINE void cotangent(const MatV & V, const MatF & F, MatC & C);
}

#ifdef IGL_HEADER_ONLY
#  include "cotangent.cpp"
#endif

#endif
