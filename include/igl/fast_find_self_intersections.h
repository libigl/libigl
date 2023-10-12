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

namespace igl
{
  /// Identify triangles where two meshes interesect 
  /// using AABBTree and tri_tri_intersection_test_3d.
  ///
  /// @param[in] V  #V by 3 list representing vertices on the first mesh
  /// @param[in] F  #F by 3 list representing triangles on the first mesh
  /// @param[out] IF #IF by 2 list of intersecting triangle pairs, so that 
  ///   F1(IF(i,0),:) intersects F2(IF(i,1),:)
  /// @param[out] EV #EV by 3 list of vertices definining intersection segments
  /// for non-coplanar intersections
  /// @param[out] EE #EE by 2 list of edges indices into rows of EV
  /// @param[out] EI #EI by 1 list of indices into rows IF indicating source of
  ///   intersection.
  ///
  /// \see copyleft::cgal::SelfIntersectMesh
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedIF,
    typename DerivedEV,
    typename DerivedEE,
    typename DerivedEI>
  IGL_INLINE bool fast_find_self_intersections(
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    const bool detect_only,
    const bool first_only,
    Eigen::PlainObjectBase<DerivedIF> & IF,
    Eigen::PlainObjectBase<DerivedEV> & EV,
    Eigen::PlainObjectBase<DerivedEE> & EE,
    Eigen::PlainObjectBase<DerivedEI> & EI);
};

#ifndef IGL_STATIC_LIBRARY
#  include "fast_find_self_intersections.cpp"
#endif

#endif
