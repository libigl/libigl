// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2021 Vladimir S. FONOV <vladimir.fonov@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/
#pragma once
#ifndef IGL_TRI_TRI_INTERSECT_H
#define IGL_TRI_TRI_INTERSECT_H

#include "igl_inline.h"
#include <Eigen/Core>

namespace igl {


// Three-dimensional Triangle-Triangle overlap test
//   if triangles are co-planar
//
// Input:
//   p1,q1,r1  - vertices of the 1st triangle (3D)
//   p2,q2,r2  - vertices of the 2nd triangle (3D)
//
// Output:
// 
//   Return true if two triangles overlap
template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2> 
IGL_INLINE bool tri_tri_overlap_test_3d(
  const Eigen::MatrixBase<DerivedP1> &  p1, 
  const Eigen::MatrixBase<DerivedQ1> &  q1, 
  const Eigen::MatrixBase<DerivedR1> &  r1, 
  const Eigen::MatrixBase<DerivedP2> &  p2, 
  const Eigen::MatrixBase<DerivedQ2> &  q2, 
  const Eigen::MatrixBase<DerivedR2> &  r2);


// Three-dimensional Triangle-Triangle Intersection Test
// additionaly computes the segment of intersection of the two triangles if it exists. 
// coplanar returns whether the triangles are coplanar, 
// source and target are the endpoints of the line segment of intersection 
//
// Input:
//   p1,q1,r1  - vertices of the 1st triangle (3D)
//   p2,q2,r2  - vertices of the 2nd triangle (3D)
//
// Output:
//   coplanar - flag if two triangles are coplanar
//   source - 1st point of intersection (if exists)
//   target - 2nd point in intersection (if exists)
//
//   Return true if two triangles intersect
template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2,
typename DerivedS,typename DerivedT>
IGL_INLINE bool tri_tri_intersection_test_3d(
    const Eigen::MatrixBase<DerivedP1> & p1, const Eigen::MatrixBase<DerivedQ1> & q1, const Eigen::MatrixBase<DerivedR1> & r1, 
    const Eigen::MatrixBase<DerivedP2> & p2, const Eigen::MatrixBase<DerivedQ2> & q2, const Eigen::MatrixBase<DerivedR2> & r2,
    bool & coplanar, 
    Eigen::MatrixBase<DerivedS> & source, 
    Eigen::MatrixBase<DerivedT> & target );



// Two dimensional Triangle-Triangle Overlap Test
// Input:
//   p1,q1,r1  - vertices of the 1st triangle (2D)
//   p2,q2,r2  - vertices of the 2nd triangle (2D)
//
// Output:
//   Return true if two triangles overlap
template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2>
IGL_INLINE bool tri_tri_overlap_test_2d(
  const Eigen::MatrixBase<DerivedP1> &p1, const Eigen::MatrixBase<DerivedQ1> &q1, const Eigen::MatrixBase<DerivedR1> &r1,
  const Eigen::MatrixBase<DerivedP2> &p2, const Eigen::MatrixBase<DerivedQ2> &q2, const Eigen::MatrixBase<DerivedR2> &r2);


};

#ifndef IGL_STATIC_LIBRARY
#  include "Guigue2003_tri_tri_intersect.cpp"
#endif

#endif // IGL_TRI_TRI_INTERSECT_H
