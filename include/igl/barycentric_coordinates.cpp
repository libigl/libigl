// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "barycentric_coordinates.h"
#include "volume.h"
#include "doublearea.h"

template <
  typename DerivedP,
  typename DerivedA,
  typename DerivedB,
  typename DerivedC,
  typename DerivedD,
  typename DerivedL>
IGL_INLINE void igl::barycentric_coordinates(
  const Eigen::PlainObjectBase<DerivedP> & P,
  const Eigen::PlainObjectBase<DerivedA> & A,
  const Eigen::PlainObjectBase<DerivedB> & B,
  const Eigen::PlainObjectBase<DerivedC> & C,
  const Eigen::PlainObjectBase<DerivedD> & D,
  Eigen::PlainObjectBase<DerivedL> & L)
{
  using namespace Eigen;
  assert(P.cols() == 3 && "query must be in 3d");
  assert(A.cols() == 3 && "corners must be in 3d");
  assert(B.cols() == 3 && "corners must be in 3d");
  assert(C.cols() == 3 && "corners must be in 3d");
  assert(D.cols() == 3 && "corners must be in 3d");
  assert(P.rows() == A.rows() && "Must have same number of queries as corners");
  assert(A.rows() == B.rows() && "Corners must be same size");
  assert(A.rows() == C.rows() && "Corners must be same size");
  assert(A.rows() == D.rows() && "Corners must be same size");
  typedef Matrix<typename DerivedL::Scalar,DerivedL::RowsAtCompileTime,1> 
    VectorXS;
  // Total volume
  VectorXS vol,LA,LB,LC,LD;
  volume(B,D,C,P,LA);
  volume(A,C,D,P,LB);
  volume(A,D,B,P,LC);
  volume(A,B,C,P,LD);
  volume(A,B,C,D,vol);
  L.resize(P.rows(),4);
  L<<LA,LB,LC,LD;
  L.array().colwise() /= vol.array();
}

template <
  typename DerivedP,
  typename DerivedA,
  typename DerivedB,
  typename DerivedC,
  typename DerivedL>
IGL_INLINE void igl::barycentric_coordinates(
  const Eigen::PlainObjectBase<DerivedP> & P,
  const Eigen::PlainObjectBase<DerivedA> & A,
  const Eigen::PlainObjectBase<DerivedB> & B,
  const Eigen::PlainObjectBase<DerivedC> & C,
  Eigen::PlainObjectBase<DerivedL> & L)
{
  using namespace Eigen;
  assert(P.cols() == 2 && "query must be in 2d");
  assert(A.cols() == 2 && "corners must be in 2d");
  assert(B.cols() == 2 && "corners must be in 2d");
  assert(C.cols() == 2 && "corners must be in 2d");
  assert(P.rows() == A.rows() && "Must have same number of queries as corners");
  assert(A.rows() == B.rows() && "Corners must be same size");
  assert(A.rows() == C.rows() && "Corners must be same size");
  typedef Matrix<typename DerivedL::Scalar,DerivedL::RowsAtCompileTime,1> 
    VectorXS;
  // Total area
  VectorXS dblA,LA,LB,LC;
  doublearea(P,B,C,LA);
  doublearea(A,P,C,LB);
  doublearea(A,B,P,LC);
  doublearea(A,B,C,dblA);
  L.resize(P.rows(),3);
  L<<LA,LB,LC;
  L.array().colwise() /= dblA.array();
}

#ifdef IGL_STATIC_LIBRARY
template void igl::barycentric_coordinates<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
