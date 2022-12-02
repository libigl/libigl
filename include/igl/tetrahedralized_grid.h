// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2020 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_TETRAHEDRALIZED_GRID_H
#define IGL_TETRAHEDRALIZED_GRID_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  enum TetrahedralizedGripType
  {
    TETRAHEDRALIZED_GRID_TYPE_5 = 0,
    TETRAHEDRALIZED_GRID_TYPE_6_ROTATIONAL = 1,
    NUM_TETRAHEDRALIZED_GRID_TYPE = 2
  };
  // Construct vertices of a regular grid, suitable for input to
  // `igl::marching_tets`
  //
  // Inputs:
  //   nx  number of grid vertices in x direction
  //   ny  number of grid vertices in y direction
  //   nz  number of grid vertices in z direction
  //   type  type of tetrahedralization of cube to use
  // Outputs:
  //   GV  nx*ny*nz by 3 list of grid vertex positions
  //   GT  (nx-1)*(ny-1)*(nz-1)*k by 4 list of tetrahedron indices into rows of
  //     V, where k is the number of tets per cube (dependent on type)
  //
  //   See also: triangulated_grid, quad_grid
  template <
    typename DerivedGV,
    typename DerivedGT>
  IGL_INLINE void tetrahedralized_grid(
    const int nx,
    const int ny,
    const int nz,
    const TetrahedralizedGripType type,
    Eigen::PlainObjectBase<DerivedGV> & GV,
    Eigen::PlainObjectBase<DerivedGT> & GT);
  //
  // Inputs:
  //   GV  nx*ny*nz by 3 list of grid vertex positions
  //   side  3-long list {nx,ny,nz} see above
  //   type  type of tetrahedralization of cube to use
  // Outputs:
  //   GT  (nx-1)*(ny-1)*(nz-1)*k by 4 list of tetrahedron indices into rows of
  //     V, where k is the number of tets per cube (dependent on type)
  //
  template <
    typename DerivedGV,
    typename Derivedside,
    typename DerivedGT>
  IGL_INLINE void tetrahedralized_grid(
    const Eigen::MatrixBase<DerivedGV> & GV,
    const Eigen::MatrixBase<Derivedside> & side,
    const TetrahedralizedGripType type,
    Eigen::PlainObjectBase<DerivedGT> & GT);
}
#ifndef IGL_STATIC_LIBRARY
#  include "tetrahedralized_grid.cpp"
#endif
#endif 

