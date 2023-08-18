// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2022 Vladimir S. FONOV <vladimir.fonov@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/
#pragma once
#ifndef FAST_FIND_SELF_INTERSECTIONS_H
#define FAST_FIND_SELF_INTERSECTIONS_H

#include "igl_inline.h"
#include <Eigen/Core>

namespace igl {

  /// Identify triangles where mesh intersects itself
  /// using AABBTree and tri_tri_intersection_test_3d
  ///
  /// @param[in] V  #V by 3 list representing vertices
  /// @param[in] F  #F by 3 list representing triangles.
  /// @param[out] intersect  #F by 1 indicator that triangle intersects anothe triangle
  /// @param[out] edges      list of pairs of intersection edges
  /// @return whether any self-interections were found
  ///
  /// \see copyleft::cgal::remesh_self_intersections
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedI,
    typename DerivedE>
  IGL_INLINE bool fast_find_self_intersections(
    const Eigen::MatrixBase<DerivedV>& V,
    const Eigen::MatrixBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedI>& intersect,
    Eigen::PlainObjectBase<DerivedE>& edges );
  /// \overload
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedI>
  IGL_INLINE bool fast_find_self_intersections(
    const Eigen::MatrixBase<DerivedV>& V,
    const Eigen::MatrixBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedI>& intersect);

};

#ifndef IGL_STATIC_LIBRARY
#  include "fast_find_self_intersections.cpp"
#endif

#endif
