// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SIGNED_DISTANCE_H
#define IGL_SIGNED_DISTANCE_H

#include "igl_inline.h"
#include "AABB.h"
#include "WindingNumberAABB.h"
#include <Eigen/Core>
#include <vector>
namespace igl
{
  enum SignedDistanceType
  {
    // Use fast pseudo-normal test [Bærentzen & Aanæs 2005]
    SIGNED_DISTANCE_TYPE_PSEUDONORMAL   = 0,
    SIGNED_DISTANCE_TYPE_WINDING_NUMBER = 1,
    SIGNED_DISTANCE_TYPE_DEFAULT        = 2,
    SIGNED_DISTANCE_TYPE_UNSIGNED       = 3,
    NUM_SIGNED_DISTANCE_TYPE            = 4
  };
  // Computes signed distance to a mesh
  //
  // Inputs:
  //   P  #P by 3 list of query point positions
  //   V  #V by 3 list of vertex positions
  //   F  #F by ss list of triangle indices, ss should be 3 unless sign_type ==
  //     SIGNED_DISTANCE_TYPE_UNSIGNED
  //   sign_type  method for computing distance _sign_ S
  // Outputs:
  //   S  #P list of smallest signed distances
  //   I  #P list of facet indices corresponding to smallest distances
  //   C  #P by 3 list of closest points
  //   N  #P by 3 list of closest normals (only set if
  //     sign_type=SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
  //
  // Known bugs: This only computes distances to triangles. So unreferenced
  // vertices and degenerate triangles are ignored.
  IGL_INLINE void signed_distance(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const SignedDistanceType sign_type,
    Eigen::VectorXd & S,
    Eigen::VectorXi & I,
    Eigen::MatrixXd & C,
    Eigen::MatrixXd & N);
  // Computes signed distance to mesh
  //
  // Inputs:
  //   tree  AABB acceleration tree (see AABB.h)
  //   F  #F by 3 list of triangle indices
  //   FN  #F by 3 list of triangle normals 
  //   VN  #V by 3 list of vertex normals (ANGLE WEIGHTING)
  //   EN  #E by 3 list of edge normals (UNIFORM WEIGHTING)
  //   EMAP  #F*3 mapping edges in F to E
  //   q  Query point
  // Returns signed distance to mesh
  //
  IGL_INLINE double signed_distance_pseudonormal(
    const AABB<Eigen::MatrixXd,3> & tree,
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & FN,
    const Eigen::MatrixXd & VN,
    const Eigen::MatrixXd & EN,
    const Eigen::VectorXi & EMAP,
    const Eigen::RowVector3d & q);
  // Outputs:
  //   s  sign
  //   sqrd  squared distance
  //   i  closest primitive
  //   c  closest point
  //   n  normal
  IGL_INLINE void signed_distance_pseudonormal(
    const AABB<Eigen::MatrixXd,3> & tree,
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & FN,
    const Eigen::MatrixXd & VN,
    const Eigen::MatrixXd & EN,
    const Eigen::VectorXi & EMAP,
    const Eigen::RowVector3d & q,
    double & s,
    double & sqrd,
    int & i,
    Eigen::RowVector3d & c,
    Eigen::RowVector3d & n);
  IGL_INLINE void signed_distance_pseudonormal(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const AABB<Eigen::MatrixXd,3> & tree,
    const Eigen::MatrixXd & FN,
    const Eigen::MatrixXd & VN,
    const Eigen::MatrixXd & EN,
    const Eigen::VectorXi & EMAP,
    Eigen::VectorXd & S,
    Eigen::VectorXi & I,
    Eigen::MatrixXd & C,
    Eigen::MatrixXd & N);

  // Inputs:
  //   tree  AABB acceleration tree (see cgal/point_mesh_squared_distance.h)
  //   hier  Winding number evaluation hierarchy
  //   q  Query point
  // Returns signed distance to mesh
  IGL_INLINE double signed_distance_winding_number(
    const AABB<Eigen::MatrixXd,3> & tree,
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const igl::WindingNumberAABB<Eigen::Vector3d> & hier,
    const Eigen::RowVector3d & q);
  // Outputs:
  //   s  sign
  //   sqrd  squared distance
  //   pp  closest point and primitve
  IGL_INLINE void signed_distance_winding_number(
    const AABB<Eigen::MatrixXd,3> & tree,
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const igl::WindingNumberAABB<Eigen::Vector3d> & hier,
    const Eigen::RowVector3d & q,
    double & s,
    double & sqrd,
    int & i,
    Eigen::RowVector3d & c);
}

#ifndef IGL_STATIC_LIBRARY
#  include "signed_distance.cpp"
#endif

#endif

