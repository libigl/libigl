// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2022 Vladimir S. FONOV <vladimir.fonov@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/
#pragma once
#ifndef FAST_FIND_MESH_INTERSECT_H
#define FAST_FIND_MESH_INTERSECT_H

#include "igl_inline.h"
#include "AABB.h"
#include <Eigen/Core>

namespace igl {
  /// Identify triangles where two meshes interesect 
  /// using AABBTree and tri_tri_intersection_test_3d.
  ///
  /// @param[in] V1  #V by 3 list representing vertices on the first mesh
  /// @param[in] F1  #F by 3 list representing triangles on the first mesh
  /// @param[in] V2  #V by 3 list representing vertices on the second mesh
  /// @param[in] F2  #F by 3 list representing triangles on the second mesh
  /// @param[out] intersect_pairs  correspondance list of intersecting triangles
  ///                    column 0 - mesh 1, column 1 - mesh2  
  /// @param[out] edges      list of pairs of intersection edges
  ///
  /// \see copyleft::cgal::intersect_other
  template <
    typename DerivedV1,
    typename DerivedF1,
    typename DerivedV2,
    typename DerivedF2,
    typename DerivedI,
    typename DerivedE>
  IGL_INLINE void fast_find_intersections(
    const Eigen::MatrixBase<DerivedV1>& V1,
    const Eigen::MatrixBase<DerivedF1>& F1,
    const Eigen::MatrixBase<DerivedV2>& V2,
    const Eigen::MatrixBase<DerivedF2>& F2,
          Eigen::PlainObjectBase<DerivedI>& intersect_pairs,
          Eigen::PlainObjectBase<DerivedE>& edges );
  /// \overload
  /// @param[in] tree - AABB tree bult from the first mesh
  template <
    typename DerivedV1,
    typename DerivedF1,
    typename DerivedV2,
    typename DerivedF2,
    typename DerivedI,
    typename DerivedE>
  IGL_INLINE void fast_find_intersections(
    const AABB<DerivedV1,3>           & tree,
    const Eigen::MatrixBase<DerivedV1>& V1,
    const Eigen::MatrixBase<DerivedF1>& F1,
    const Eigen::MatrixBase<DerivedV2>& V2,
    const Eigen::MatrixBase<DerivedF2>& F2,
          Eigen::PlainObjectBase<DerivedI>& intersect_pairs,
          Eigen::PlainObjectBase<DerivedE>& edges );
};

#ifndef IGL_STATIC_LIBRARY
#  include "fast_find_intersections.cpp"
#endif

#endif
