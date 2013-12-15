// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "polar_svd.h"
#include <Eigen/SVD>

// Adapted from Olga's CGAL mentee's ARAP code
template <
  typename DerivedA,
  typename DerivedR,
  typename DerivedT>
IGL_INLINE void igl::polar_svd(
  const Eigen::PlainObjectBase<DerivedA> & A,
  Eigen::PlainObjectBase<DerivedR> & R,
  Eigen::PlainObjectBase<DerivedT> & T)
{
  Eigen::PlainObjectBase<DerivedA> U;
  Eigen::PlainObjectBase<DerivedA> V;
  Eigen::Matrix<typename DerivedA::Scalar,DerivedA::RowsAtCompileTime,1> S;
  return igl::polar_svd(A,R,T,U,S,V);
}

template <
  typename DerivedA,
  typename DerivedR,
  typename DerivedT,
  typename DerivedU,
  typename DerivedS,
  typename DerivedV>
IGL_INLINE void igl::polar_svd(
  const Eigen::PlainObjectBase<DerivedA> & A,
  Eigen::PlainObjectBase<DerivedR> & R,
  Eigen::PlainObjectBase<DerivedT> & T,
  Eigen::PlainObjectBase<DerivedU> & U,
  Eigen::PlainObjectBase<DerivedS> & S,
  Eigen::PlainObjectBase<DerivedV> & V)
{
  Eigen::JacobiSVD<DerivedA> svd;
  svd.compute(A, Eigen::ComputeFullU | Eigen::ComputeFullV );
  U = svd.matrixU();
  V = svd.matrixV();
  S = svd.singularValues();
  R = U*V.transpose();
  T = V*S.asDiagonal()*V.adjoint();
}

#ifndef IGL_HEADER_ONLY
// Explicit template instanciation
template void igl::polar_svd<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::polar_svd<Eigen::Matrix<double, 2, 2, 0, 2, 2>, Eigen::Matrix<double, 2, 2, 0, 2, 2>, Eigen::Matrix<double, 2, 2, 0, 2, 2> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 2, 2, 0, 2, 2> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 2, 2, 0, 2, 2> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, 2, 2, 0, 2, 2> >&);
template void igl::polar_svd<Eigen::Matrix<double, 3, 3, 0, 3, 3>, Eigen::Matrix<double, 3, 3, 0, 3, 3>, Eigen::Matrix<double, 3, 3, 0, 3, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 3, 0, 3, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 3, 0, 3, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 3, 0, 3, 3> >&);
template void igl::polar_svd<Eigen::Matrix<float, 2, 2, 0, 2, 2>, Eigen::Matrix<float, 2, 2, 0, 2, 2>, Eigen::Matrix<float, 2, 2, 0, 2, 2> >(Eigen::PlainObjectBase<Eigen::Matrix<float, 2, 2, 0, 2, 2> > const&, Eigen::PlainObjectBase<Eigen::Matrix<float, 2, 2, 0, 2, 2> >&, Eigen::PlainObjectBase<Eigen::Matrix<float, 2, 2, 0, 2, 2> >&);
#endif
