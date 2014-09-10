// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SIGNED_DISTANCE_H
#define IGL_SIGNED_DISTANCE_H
#include <igl/igl_inline.h>
#include <igl/WindingNumberAABB.h>
#include <Eigen/Core>
#include <vector>
#include "CGAL_includes.hpp"
namespace igl
{
  enum SignedDistanceType
  {
    SIGNED_DISTANCE_TYPE_PSEUDONORMAL   = 0,
    SIGNED_DISTANCE_TYPE_WINDING_NUMBER = 1,
    SIGNED_DISTANCE_TYPE_DEFAULT        = 2,
    NUM_SIGNED_DISTANCE_TYPE            = 3
  };
  // Computes signed distance to mesh
  //
  // Templates:
  //   Kernal  CGAL computation and construction kernel (e.g.
  //     CGAL::Simple_cartesian<double>)
  // Inputs:
  //   tree  AABB acceleration tree (see point_mesh_squared_distance.h)
  //   T  #F list of CGAL triangles (see point_mesh_squared_distance.h)
  //   F  #F by 3 list of triangle indices
  //   FN  #F by 3 list of triangle normals 
  //   VN  #V by 3 list of vertex normals (ANGLE WEIGHTING)
  //   EN  #E by 3 list of edge normals (UNIFORM WEIGHTING)
  //   EMAP  #F*3 mapping edges in F to E
  //   q  Query point
  // Returns signed distance to mesh
  //
  template <typename Kernel>
  IGL_INLINE typename Kernel::FT signed_distance_pseudonormal(
    const CGAL::AABB_tree<
      CGAL::AABB_traits<Kernel, 
        CGAL::AABB_triangle_primitive<Kernel, 
          typename std::vector<CGAL::Triangle_3<Kernel> >::iterator
        >
      >
    > & tree,
    const std::vector<CGAL::Triangle_3<Kernel> > & T,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & FN,
    const Eigen::MatrixXd & VN,
    const Eigen::MatrixXd & EN,
    const Eigen::VectorXi & EMAP,
    const typename Kernel::Point_3 & q);
  // Inputs:
  //   tree  AABB acceleration tree (see point_mesh_squared_distance.h)
  //   hier  Winding number evaluation hierarchy
  //   q  Query point
  // Returns signed distance to mesh
  template <typename Kernel>
  IGL_INLINE typename Kernel::FT signed_distance_winding_number(
    const CGAL::AABB_tree<
      CGAL::AABB_traits<Kernel, 
        CGAL::AABB_triangle_primitive<Kernel, 
          typename std::vector<CGAL::Triangle_3<Kernel> >::iterator
        >
      >
    > & tree,
    const igl::WindingNumberAABB<Eigen::Vector3d> & hier,
    const typename Kernel::Point_3 & q);
}

#ifndef IGL_STATIC_LIBRARY
#  include "signed_distance.cpp"
#endif

#endif

