// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Daniele Panozzo <daniele.panozzo@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "line_mesh_intersection.h"

// For error printing
#include <cstdio>
#include <vector>

#include <igl/per_vertex_normals.h>
#include <igl/embree/EmbreeIntersector.h>

//template <typename ScalarMatrix, typename IndexMatrix>
//IGL_INLINE ScalarMatrix igl::line_mesh_intersection(
//   const ScalarMatrix & V_source,
//   const IndexMatrix  & F_source,
//   const ScalarMatrix & V_target,
//   const IndexMatrix  & F_target
//)
//{
//  // Compute normals for the tri
//  Eigen::MatrixXd ray_dir;
//  igl::per_vertex_normals(V_source, F_source, ray_dir);
//
//  return line_mesh_intersection(V_source,ray_dir,V_target,F_target);
//}

template <typename ScalarMatrix, typename IndexMatrix>
IGL_INLINE ScalarMatrix igl::line_mesh_intersection
(
 const ScalarMatrix & V_source,
 const ScalarMatrix  & N_source,
 const ScalarMatrix & V_target,
 const IndexMatrix  & F_target
)
{

  double tol = 0.00001;
  
  Eigen::MatrixXd ray_pos = V_source;
  Eigen::MatrixXd ray_dir = N_source;

  // Allocate matrix for the result
  ScalarMatrix R;
  R.resize(V_source.rows(), 3);
  
  // Initialize embree
  igl::EmbreeIntersector embree;
  embree.init(V_target.template cast<float>(),F_target.template cast<int>());

  // Shoot rays from the source to the target
  for (unsigned i=0; i<ray_pos.rows(); ++i)
  {
    igl::Hit A,B;
    
    // Shoot ray A
    Eigen::RowVector3d A_pos = ray_pos.row(i) + tol * ray_dir.row(i);
    Eigen::RowVector3d A_dir = -ray_dir.row(i);
    
    bool A_hit = embree.intersectRay(A_pos.cast<float>(), A_dir.cast<float>(),A);
    
    Eigen::RowVector3d B_pos = ray_pos.row(i) - tol * ray_dir.row(i);
    Eigen::RowVector3d B_dir = ray_dir.row(i);
    
    bool B_hit = embree.intersectRay(B_pos.cast<float>(), B_dir.cast<float>(),B);
    
    
    int choice = -1;
    
    if (A_hit && ! B_hit)
      choice = 0;
    else if (!A_hit && B_hit)
      choice = 1;
    else if (A_hit && B_hit)
      choice = A.t > B.t;
    
    Eigen::RowVector3d temp;

    if (choice == -1)
      temp << -1, 0, 0;
    else if (choice == 0)
      temp << A.id, A.u, A.v;
    else if (choice == 1)
      temp << B.id, B.u, B.v;

    R.row(i) = temp;
    
  }

  return R;

}

