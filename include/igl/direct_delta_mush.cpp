// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2020 Xiangyu Kong <xiangyu.kong@mail.utoronto.ca>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "direct_delta_mush.h"

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedC,
  typename DerivedE,
  typename DerivedT,
  typename DerivedU>
IGL_INLINE void direct_delta_mush(
  const Eigen::MatrixBase<DerivedV> &V,
  const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedC> &C,
  const Eigen::MatrixBase<DerivedE> &E,
  const std::vector<
    Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>> &T,
  Eigen::PlainObjectBase<DerivedU> &U)
{
  return;
}

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedC,
  typename DerivedE,
  typename DerivedOmega>
IGL_INLINE void precomputation(
  const Eigen::MatrixBase<DerivedV> &V,
  const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedC> &C,
  const Eigen::MatrixBase<DerivedE> &E,
  const int p,
  const typename DerivedV::Scalar lambda,
  const typename DerivedV::Scalar kappa,
  Eigen::PlainObjectBase<DerivedOmega> &Omega)
{
  return;
}

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedC,
  typename DerivedE,
  typename DerivedT,
  typename DerivedOmega,
  typename DerivedU>
IGL_INLINE void pose_evaluation(
  const std::vector<
    Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d>> &T,
  const Eigen::MatrixBase<DerivedOmega> &Omega,
  Eigen::PlainObjectBase<DerivedU> &U)
{
  return;
}