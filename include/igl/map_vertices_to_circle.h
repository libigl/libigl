// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Stefan Brugger <stefanbrugger@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_MAPVERTICESTOCIRCLE_H
#define IGL_MAPVERTICESTOCIRCLE_H
#include <igl/igl_inline.h>

#include <Eigen/Dense>
#include <vector>

namespace igl
{

  // Map the vertices whose indices are in b on the unit circle.
  //
  // Inputs:
  //   V  #V by dim list of mesh vertex positions
  //   F  #V by dim list of mesh faces
  //   b  #W list of vertex ids
  // Outputs:
  //   UV   #W by 2 list of 2D position on the unit circle for the vertices in b
  IGL_INLINE void map_vertices_to_circle(
  	const Eigen::MatrixXd& V,
  	const Eigen::MatrixXi& F,
    const Eigen::VectorXi& b,
  	Eigen::MatrixXd& UV);
}

#ifdef IGL_HEADER_ONLY
#  include "map_vertices_to_circle.cpp"
#endif

#endif
