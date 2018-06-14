// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Gavin Barill <gavinpcb@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/

#ifndef IGL_BUILD_OCTREE
#define IGL_BUILD_OCTREE
#include <Eigen/Core>
#include "igl_inline.h"

namespace igl
{
  // Given a set of 3D points P, generate data structures for a pointerless
  // octree. Each cell stores its points, children, center location and width.
  // Our octree is not dense. We use the following rule: if the current cell
  // has any number of points, it will have all 8 children. A leaf cell will
  // have -1's as its list of child indices.
  //
  // We use a binary numbering of children. Treating the parent cell's center
  // as the origin, we number the octants in the following manner:
  // The first bit is 1 iff the octant's x coordinate is positive
  // The second bit is 1 iff the octant's y coordinate is positive
  // The third bit is 1 iff the octant's z coordinate is positive
  //
  // For example, the octant with negative x, positive y, positive z is:
  // 110 binary = 6 decimal
  //
  // Inputs:
  //   P  #P by 3 list of point locations
  //
  // Outputs:
  //   point_indices  a vector of vectors, where the ith entry is a vector of
  //                  the indices into P that are the ith octree cell's points
  //   children       a vector of vectors, where the ith entry is a vector of
  //                  the ith octree cell's of octree children
  //   centers        a vector where the ith entry is a 3d row vector
  //                  representing the position of the ith cell's center
  //   widths          a vector where the ith entry is the width of the ith
  //                  octree cell
  //
  template <typename DerivedP, typename IndexType, typename CentersType,
    typename WidthsType>
  IGL_INLINE void build_octree(const Eigen::MatrixBase<DerivedP>& P,
    std::vector<std::vector<IndexType> > & point_indices,
    std::vector<Eigen::Matrix<IndexType,8,1>,
      Eigen::aligned_allocator<Eigen::Matrix<IndexType,8,1> > > & children,
    std::vector<Eigen::Matrix<CentersType,1,3>,
      Eigen::aligned_allocator<Eigen::Matrix<CentersType,1,3> > > & centers,
    std::vector<WidthsType> & widths);
}

#ifndef IGL_STATIC_LIBRARY
#  include "build_octree.cpp"
#endif

#endif

