// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2020 Xiangyu Kong <xiangyu.kong@mail.utoronto.ca>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "direct_delta_mush.h"
#include "cotmatrix.h"
#include "diag.h"


// ===== DEBUG: START
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/LU>
#include <iostream>

using namespace std;
// ===== DEBUG: END

// TODOs:
// 1. All `Scalar`s are temporarily substituted as `DerivedV`

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedC,
  typename DerivedE,
  typename DerivedW,
  typename DerivedT,
  typename DerivedTAlloc,
  typename DerivedU>
IGL_INLINE void igl::direct_delta_mush(
  const Eigen::MatrixBase<DerivedV> &V,
  const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedC> &C,
  const Eigen::MatrixBase<DerivedE> &E,
  const Eigen::MatrixBase<DerivedW> &W,
  const std::vector<DerivedT, DerivedTAlloc> &T,
  Eigen::PlainObjectBase<DerivedU> &U)
{
  cout << "START DDM" << endl;
  cout << "END DDM" << endl;
}

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedC,
  typename DerivedE,
  typename DerivedW,
  typename DerivedOmega>
IGL_INLINE void igl::direct_delta_mush_precomputation(
  const Eigen::MatrixBase<DerivedV> &V,
  const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedC> &C,
  const Eigen::MatrixBase<DerivedE> &E,
  const Eigen::MatrixBase<DerivedW> &W,
  const int p,
  const typename DerivedV::Scalar lambda,
  const typename DerivedV::Scalar kappa,
  Eigen::PlainObjectBase<DerivedOmega> &Omega)
{
  cout << "START DDM Precomputation" << endl;
  cout << "END DDM Precomputation" << endl;
}

template <
  typename DerivedT,
  typename DerivedTAlloc,
  typename DerivedOmega,
  typename DerivedU>
IGL_INLINE void igl::direct_delta_mush_pose_evaluation(
  const std::vector<DerivedT, DerivedTAlloc> &T,
  const Eigen::MatrixBase<DerivedOmega> &Omega,
  Eigen::PlainObjectBase<DerivedU> &U)
{
  cout << "START DDM Pose Eval" << endl;
  cout << "END DDM Pose Eval" << endl;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::direct_delta_mush<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Transform<double, 3, 2, 0>, Eigen::aligned_allocator<Eigen::Transform<double, 3, 2, 0> >, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, std::__1::vector<Eigen::Transform<double, 3, 2, 0>, Eigen::aligned_allocator<Eigen::Transform<double, 3, 2, 0> > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::direct_delta_mush_precomputation<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, int, Eigen::Matrix<double, -1, -1, 0, -1, -1>::Scalar, Eigen::Matrix<double, -1, -1, 0, -1, -1>::Scalar, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif