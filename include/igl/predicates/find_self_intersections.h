// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2024 Alec Jacobson alecjacobson@gmail.com
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/
#ifndef IGL_PREDICATES_FIND_SELF_INTERSECTIONS_H
#define IGL_PREDICATES_FIND_SELF_INTERSECTIONS_H

#include "../igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  namespace predicates
  {
    /// Identify triangle-triangle interesections within the same mesh using
    /// AABBTree and igl::predicates::triangle_triangle_intersect.
    ///
    /// @param[in] V  #V by 3 list representing vertices on the first mesh
    /// @param[in] F  #F by 3 list representing triangles on the first mesh
    /// @param[out] IF #IF by 2 list of intersecting triangle pairs, so that 
    ///   F1(IF(i,0),:) intersects F2(IF(i,1),:)
    ///
    /// \see copyleft::cgal::SelfIntersectMesh
    template <
      typename DerivedV,
      typename DerivedF,
      typename DerivedIF>
    IGL_INLINE bool find_self_intersections(
      const Eigen::MatrixBase<DerivedV> & V,
      const Eigen::MatrixBase<DerivedF> & F,
      const bool first_only,
      Eigen::PlainObjectBase<DerivedIF> & IF);
  }
}

#ifndef IGL_STATIC_LIBRARY
#  include "find_self_intersections.cpp"
#endif

#endif

