// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "face_areas.h"
#include "edge_lengths.h"
#include "doublearea.h"

template <typename DerivedV, typename DerivedT, typename DerivedA>
IGL_INLINE void igl::face_areas(
  const Eigen::MatrixBase<DerivedV>& V,
  const Eigen::MatrixBase<DerivedT>& T,
  Eigen::PlainObjectBase<DerivedA>& A)
{
  assert(T.cols() == 4);
  DerivedA L;
  edge_lengths(V,T,L);
  return face_areas(L,A);
}

template <typename DerivedL, typename DerivedA>
IGL_INLINE void igl::face_areas(
  const Eigen::MatrixBase<DerivedL>& L,
  Eigen::PlainObjectBase<DerivedA>& A)
{
  return face_areas(
    L,std::numeric_limits<typename DerivedL::Scalar>::quiet_NaN(),A);
}

template <typename DerivedL, typename DerivedA>
IGL_INLINE void igl::face_areas(
  const Eigen::MatrixBase<DerivedL>& L,
  const typename DerivedL::Scalar doublearea_nan_replacement,
  Eigen::PlainObjectBase<DerivedA>& A)
{
  assert(L.cols() == 6);
  const int m = L.rows();
  // (unsigned) face Areas (opposite vertices: 1 2 3 4)
  Eigen::Matrix<typename DerivedA::Scalar ,Eigen::Dynamic,1>
    A0(m,1), A1(m,1), A2(m,1), A3(m,1);
  Eigen::Matrix<typename DerivedA::Scalar ,Eigen::Dynamic,3>
    L0(m,3), L1(m,3), L2(m,3), L3(m,3);
  L0<<L.col(1),L.col(2),L.col(3);
  L1<<L.col(0),L.col(2),L.col(4);
  L2<<L.col(0),L.col(1),L.col(5);
  L3<<L.col(3),L.col(4),L.col(5);
  doublearea(L0,doublearea_nan_replacement,A0);
  doublearea(L1,doublearea_nan_replacement,A1);
  doublearea(L2,doublearea_nan_replacement,A2);
  doublearea(L3,doublearea_nan_replacement,A3);
  A.resize(m,4);
  A.col(0) = 0.5*A0;
  A.col(1) = 0.5*A1;
  A.col(2) = 0.5*A2;
  A.col(3) = 0.5*A3;
}
#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::face_areas<Eigen::Matrix<double, -1, 6, 0, -1, 6>, Eigen::Matrix<double, -1, 4, 0, -1, 4> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, 6, 0, -1, 6> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 4, 0, -1, 4> >&);
#endif
