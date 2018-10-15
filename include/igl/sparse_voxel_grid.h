// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Francis Williams <francis@fwilliams.info>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SPARSE_VOXEL_GRID_H
#define IGL_SPARSE_VOXEL_GRID_H

#include "igl_inline.h"

#include <Eigen/Core>

namespace igl {

  // sparse_voxel_grid( p0, scalarFunc, eps, CV, CS, CI )
  //
  // Given a point, p0, on an isosurface, construct a shell of epsilon sized cubes surrounding the surface.
  // These cubes can be used as the input to marching cubes.
  //
  // Input:
  //  p0  A 3D point on the isosurface surface defined by scalarFunc(x) = 0
  //  scalarFunc  A scalar function from R^3 to R -- points which map to 0 lie
  //              on the surface, points which are negative lie inside the surface,
  //              and points which are positive lie outside the surface
  //  eps  The edge length of the cubes surrounding the surface
  //  expected_number_of_cubes  This pre-allocates internal data structures to speed things up
  // Output:
  //   CS  #cube-vertices by 1 list of scalar values at the cube vertices
  //   CV  #cube-vertices by 3 list of cube vertex positions
  //   CI  #number of cubes by 8 list of indexes into CS and CV. Each row represents a cube
  //
  template <typename DerivedS, typename DerivedP0, typename DerivedV, typename DerivedI>
  IGL_INLINE void sparse_voxel_grid(const Eigen::MatrixBase<DerivedP0>& p0,
                                    const std::function<typename DerivedS::Scalar(const DerivedP0&)>& scalarFunc,
                                    const double eps,
                                    Eigen::PlainObjectBase<DerivedS>& CS,
                                    Eigen::PlainObjectBase<DerivedV>& CV,
                                    Eigen::PlainObjectBase<DerivedI>& CI,
                                    int expected_number_of_cubes=1024);

}
#ifndef IGL_STATIC_LIBRARY
#    include "sparse_voxel_grid.cpp"
#endif

#endif // IGL_SPARSE_VOXEL_GRID_H
