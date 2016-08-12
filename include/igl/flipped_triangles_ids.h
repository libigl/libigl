// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Michael Rabinovich <michaelrabinovich27@gmail.com@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_FLIPPED_TRIANGLES_IDS_H
#define IGL_FLIPPED_TRIANGLES_IDS_H
#include "igl_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
namespace igl
{
  // Finds the ids of the flipped triangles of the mesh V,F in the UV mapping uv
  // Templates:
  //   Scalar  should be a floating point number type
  //   Index   should be an integer type
  // Inputs:
  //   V       #V by 2 list of mesh vertex positions
  //   F       #F by dim list of mesh faces (must be triangles)
  // Outputs:
  //   A vector containing the ids of the flipped triangles
  template <typename Scalar, typename Index>
  IGL_INLINE Eigen::VectorXi flipped_triangles_ids(
    const Eigen::PlainObjectBase<Scalar> & V,
    const Eigen::PlainObjectBase<Index> & F
  );

}

#ifndef IGL_STATIC_LIBRARY
#  include "flipped_triangles_ids.cpp"
#endif

#endif
