// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Gavin Barill <gavinpcb@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/

#ifndef IGL_KNN_OCTREE
#define IGL_KNN_OCTREE
#include "igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  // Given a 3D set of points P, an whole number k, and an octree
  // find the indicies of the k nearest neighbors for each point in P.
  // Note that each point is its own neighbor.
  //
  // The octree data structures used in this function are intended to be the
  // same ones output from igl::build_octree
  //
  // Inputs:
  //   P  #P by 3 list of point locations
  //   k  number of neighbors to find
  //   point_indices  a vector of vectors, where the ith entry is a vector of
  //                  the indices into P that are the ith octree cell's points
  //   children       a vector of vectors, where the ith entry is a vector of
  //                  the ith octree cell's of octree children
  //   centers        a vector where the ith entry is a 3d row vector
  //                  representing the position of the ith cell's center
  //   widths         a vector where the ith entry is the width of the ith
  //                  octree cell
  // Outputs:
  //   I  #P by k list of k-nearest-neighbor indices into P
  template <typename DerivedP, typename KType, typename IndexType,
    typename CentersType, typename WidthsType, typename DerivedI>
  IGL_INLINE void knn_octree(const Eigen::MatrixBase<DerivedP>& P,
    const KType & k,
    const std::vector<std::vector<IndexType> > & point_indices,
    const std::vector<Eigen::Matrix<IndexType,8,1>, Eigen::aligned_allocator<Eigen::Matrix<IndexType,8,1> > > & children,
    const std::vector<Eigen::Matrix<CentersType,1,3>, Eigen::aligned_allocator<Eigen::Matrix<CentersType,1,3> > > & centers,
    const std::vector<WidthsType> & widths,
    Eigen::PlainObjectBase<DerivedI> & I);
}
#ifndef IGL_STATIC_LIBRARY
#  include "knn_octree.cpp"
#endif
#endif

