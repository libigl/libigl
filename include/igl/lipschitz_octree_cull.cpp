// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2025 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "lipschitz_octree_cull.h"
#include "unique_sparse_voxel_corners.h"
#include "find.h"
#include "matlab_format.h"
#include "parallel_for.h"
#include <cassert>
#include <iostream>


template <
  typename Derivedorigin,
  typename Func,
  typename Derivedijk,
  typename Derivedijk_maybe
    >
IGL_INLINE void igl::lipschitz_octree_cull(
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

  Eigen::VectorXi I;
  Eigen::Matrix<int,Eigen::Dynamic,8,Eigen::RowMajor> J;
  MatrixiX3R unique_ijk;
  MatrixSX3R unique_corner_positions;
  igl::unique_sparse_voxel_corners(origin,h0,depth,ijk,unique_ijk,I,J,unique_corner_positions);

  Eigen::Matrix<Scalar,Eigen::Dynamic,1> uU(unique_corner_positions.rows());
  //for(int u = 0;u<unique_corner_positions.rows();u++)
  igl::parallel_for(
    unique_corner_positions.rows(),
    [&](const int u)
    {
      // evaluate the function at the corner
      const RowVectorS3 corner = unique_corner_positions.row(u);
      uU(u) = udf(corner);
      assert(uU(u) >= 0 && "udf must be non-negative for lipschitz_octree_cull");
    },
    1000);

  MatrixSX8R U(ijk.rows(),8);
  // Pull this out as a function that identifies unique corners and then batches
  // the udf calls.
  // Should use parallel_for
  for(int c = 0;c<ijk.rows();c++)
  {
    for(int i = 0;i<8;i++)
    {
      U(c,i) = uU(J(c,i));
    }
  }
  // It's a bit silly to expose an overload for the _replicated_ U values. The
  // subsequent test is done elementwize so that could have been done before
  // mapping back to replicates.
  return lipschitz_octree_cull(h, U, ijk, ijk_maybe);
}

template <
  typename DerivedU,
  typename Derivedijk,
  typename Derivedijk_maybe
    >
IGL_INLINE void igl::lipschitz_octree_cull(
  const typename DerivedU::Scalar h,
  const Eigen::MatrixBase<DerivedU> & U,
  const Eigen::MatrixBase<Derivedijk> & ijk,
  Eigen::PlainObjectBase<Derivedijk_maybe> & ijk_maybe)
{
  assert(U.rows() == ijk.rows() &&
    "U and ijk must have the same number of rows");
  assert(U.cols() == 8 &&
    "U must have 8 columns for the 8 corners of each octree cell");
  assert(ijk.cols() == 3 &&
    "ijk must have 3 columns for x, y, z indices");
  using MatrixiX3R = Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor>;
  // sufficient_1 = any(U > h*sqrt(3),2);
  using ArraybX = Eigen::Array<bool,Eigen::Dynamic,1>;
  const ArraybX empty = (U.array() > h * std::sqrt(3)).rowwise().any().eval();
  // maybe = !empty
  const ArraybX maybe = (!empty).eval();
  //const auto maybe = (U.array() < h * std::sqrt(3)).rowwise().all().eval();
  ijk_maybe = ijk(igl::find(maybe), Eigen::all);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::lipschitz_octree_cull<Eigen::Matrix<double, 1, 3, 1, 1, 3>, std::function<double (Eigen::Matrix<double, 1, 3, 1, 1, 3> const&)>, Eigen::Matrix<int, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 3, 1, -1, 3>>(Eigen::MatrixBase<Eigen::Matrix<double, 1, 3, 1, 1, 3>> const&, Eigen::Matrix<double, 1, 3, 1, 1, 3>::Scalar, int, std::function<double (Eigen::Matrix<double, 1, 3, 1, 1, 3> const&)> const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3>> const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 1, -1, 3>>&);
#endif

