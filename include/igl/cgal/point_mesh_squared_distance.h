// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_POINT_MESH_SQUARED_DISTANCE_H
#define IGL_POINT_MESH_SQUARED_DISTANCE_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include "CGAL_includes.hpp"
namespace igl
{
  // Compute distances from a set of points P to a triangle mesh (V,F)
  //
  // Templates:
  //   Kernal  CGAL computation and construction kernel (e.g.
  //     CGAL::Simple_cartesian<double>)
  // Inputs:
  //   P  #P by 3 list of query point positions
  //   V  #V by 3 list of vertex positions
  //   F  #F by 3 list of triangle indices
  // Outputs:
  //   sqrD  #P list of smallest squared distances
  //   I  #P list of facet indices corresponding to smallest distances
  //   C  #P by 3 list of closest points
  //
  // Known bugs: This only computes distances to triangles. So unreferenced
  // vertices are ignored.
  template <typename Kernel>
  IGL_INLINE void point_mesh_squared_distance(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::VectorXd & sqrD,
    Eigen::VectorXi & I,
    Eigen::MatrixXd & C);
}
#ifndef IGL_STATIC_LIBRARY
#  include "point_mesh_squared_distance.cpp"
#endif

#endif
