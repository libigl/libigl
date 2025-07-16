// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2025 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "lipschitz_octree_prune.h"
#include "unique_sparse_voxel_corners.h"
#include "find.h"
#include "matlab_format.h"
#include "parallel_for.h"
#include <cassert>
#include <iostream>


template <
  bool batched,
  typename Derivedorigin,
  typename Func,
  typename Derivedijk,
  typename Derivedijk_maybe
    >
IGL_INLINE void igl::lipschitz_octree_prune(
  const Eigen::MatrixBase<Derivedorigin> & origin,
  const typename Derivedorigin::Scalar h0,
  const int depth,
  const Func & udf,
  const Eigen::MatrixBase<Derivedijk> & ijk,
  Eigen::PlainObjectBase<Derivedijk_maybe> & ijk_maybe)
{
  // static assert to ensure that Derivedorigin is a vector and the
  // non-singleton dimension is 3 or Eigen::Dynamic
  static_assert(
    (Derivedorigin::RowsAtCompileTime == 1 && (
      Derivedorigin::ColsAtCompileTime == 3 ||
      Derivedorigin::ColsAtCompileTime == Eigen::Dynamic)) ||
    (Derivedorigin::ColsAtCompileTime == 1 && (
      Derivedorigin::RowsAtCompileTime == 3 ||
      Derivedorigin::RowsAtCompileTime == Eigen::Dynamic)),
    "Derivedorigin must be a vector with 3 or Eigen::Dynamic dimensions");
  // dynamic assert that the origin is a 3D vector
  assert((origin.rows() == 3 || origin.cols() == 3) && origin.size() == 3 &&
    "origin must be a 3D vector");

  using Scalar = typename Derivedorigin::Scalar;
  using RowVectorS3 = Eigen::Matrix<Scalar,1,3>;
  using MatrixSX8R = Eigen::Matrix<typename Derivedorigin::Scalar,Eigen::Dynamic,8,Eigen::RowMajor>;
  using MatrixSX3R = Eigen::Matrix<typename Derivedorigin::Scalar,Eigen::Dynamic,3,Eigen::RowMajor>;
  using MatrixiX3R = Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor>;
  // h0 is already the h at this depth.
  const Scalar h = h0 / (1 << depth);

  Eigen::Matrix<int,Eigen::Dynamic,8,Eigen::RowMajor> J;
  MatrixiX3R unique_ijk;
  MatrixSX3R unique_corner_positions;
  igl::unique_sparse_voxel_corners(origin,h0,depth,ijk,unique_ijk,J,unique_corner_positions);

  // Effectively a batched call to udf
  Eigen::Array<bool,Eigen::Dynamic,1> big(unique_corner_positions.rows()); 
  // Requires C++17
  if constexpr (batched)
  {
    Eigen::Matrix<Scalar,Eigen::Dynamic,1> u = udf(unique_corner_positions);
    for(int i = 0;i<unique_corner_positions.rows();i++)
    {
      assert(u(i) >= 0 && "udf must be non-negative for lipschitz_octree_prune");
      big(i) = (u(i) > h * std::sqrt(3));
    }
  }else
  {
    //for(int u = 0;u<unique_corner_positions.rows();u++)
    igl::parallel_for(
      unique_corner_positions.rows(),
      [&](const int i)
      {
        // evaluate the function at the corner
        const RowVectorS3 corner = unique_corner_positions.row(i);
        const Scalar u = udf(corner);
        assert(u >= 0 && "udf must be non-negative for lipschitz_octree_prune");
        big(i) = (u > h * std::sqrt(3));
      },
      1000);
  }

  ijk_maybe.resize(ijk.rows(),3);
  int k = 0;
  for(int c = 0;c<ijk.rows();c++)
  {
    bool empty = false;
    for(int i = 0;i<8;i++)
    {
      empty = empty || big(J(c,i));
    }
    bool maybe = !empty;
    if(maybe)
    {
      ijk_maybe.row(k++) = ijk.row(c);
    }
  }
  ijk_maybe.conservativeResize(k,3);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::lipschitz_octree_prune<false,Eigen::Matrix<double, 1, 3, 1, 1, 3>, std::function<double (Eigen::Matrix<double, 1, 3, 1, 1, 3> const&)>, Eigen::Matrix<int, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 3, 1, -1, 3>>(Eigen::MatrixBase<Eigen::Matrix<double, 1, 3, 1, 1, 3>> const&, Eigen::Matrix<double, 1, 3, 1, 1, 3>::Scalar, int, std::function<double (Eigen::Matrix<double, 1, 3, 1, 1, 3> const&)> const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3>> const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 1, -1, 3>>&);
template void igl::lipschitz_octree_prune<true, Eigen::Matrix<double, 1, 3, 1, 1, 3>, std::function<Eigen::Matrix<double, -1, 1, 0, -1, 1> (Eigen::Matrix<double, -1, 3, 1, -1, 3> const&)>, Eigen::Matrix<int, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 3, 1, -1, 3>>(Eigen::MatrixBase<Eigen::Matrix<double, 1, 3, 1, 1, 3>> const&, Eigen::Matrix<double, 1, 3, 1, 1, 3>::Scalar, int, std::function<Eigen::Matrix<double, -1, 1, 0, -1, 1> (Eigen::Matrix<double, -1, 3, 1, -1, 3> const&)> const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3>> const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 1, -1, 3>>&);
#endif

