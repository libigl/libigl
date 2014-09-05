// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "internal_angles.h"
#include "edge_lengths.h"

template <typename DerivedV, typename DerivedF, typename DerivedK>
IGL_INLINE void igl::internal_angles(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedF>& F,
  Eigen::PlainObjectBase<DerivedK> & K)
{
  using namespace Eigen;
  // Edge lengths
  Matrix<
    typename DerivedV::Scalar,
    DerivedF::RowsAtCompileTime,
    DerivedF::ColsAtCompileTime> L;
  edge_lengths(V,F,L);

  assert(F.cols() == 3 && "F should contain triangles");
  return internal_angles(L,K);
}

template <typename DerivedL, typename DerivedK>
IGL_INLINE void igl::internal_angles(
  const Eigen::PlainObjectBase<DerivedL>& L,
  Eigen::PlainObjectBase<DerivedK> & K)
{
  assert(L.cols() == 3 && "Edge-lengths should come from triangles");
  K.resize(L.rows(),L.cols());
  for(int d = 0;d<3;d++)
  {
    const auto & s1 = L.col(d).array();
    const auto & s2 = L.col((d+1)%3).array();
    const auto & s3 = L.col((d+2)%3).array();
    K.col(d) = ((s3.square() + s2.square() - s1.square())/(2.*s3*s2)).acos();
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::internal_angles<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::internal_angles<Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<unsigned int, -1, -1, 1, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<unsigned int, -1, -1, 1, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
