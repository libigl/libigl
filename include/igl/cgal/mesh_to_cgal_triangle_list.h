// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_MESH_TO_CGAL_TRIANGLE_LIST_H
#define IGL_MESH_TO_CGAL_TRIANGLE_LIST_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include "CGAL_includes.hpp"
namespace igl
{
  // Convert a mesh (V,F) to a list of CGAL triangles
  //
  // Templates:
  //   Kernal  CGAL computation and construction kernel (e.g.
  //     CGAL::Exact_predicates_exact_constructions_kernel)
  // Inputs:
  //   V  #V by 3 list of vertex positions
  //   F  #F by 3 list of triangle indices
  // Outputs:
  //   T  #F list of CGAL triangles
  template <typename Kernel>
  IGL_INLINE void mesh_to_cgal_triangle_list(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    std::vector<CGAL::Triangle_3<Kernel> > & T);
}
#ifndef IGL_STATIC_LIBRARY
#  include "mesh_to_cgal_triangle_list.cpp"
#endif

#endif
