// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2018 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "accumarray.h"
#include <cassert>

template <
  typename DerivedS,
  typename DerivedV,
  typename DerivedA
  >
void igl::accumarray(
  const Eigen::MatrixBase<DerivedS> & S,
  const Eigen::MatrixBase<DerivedV> & V,
  Eigen::PlainObjectBase<DerivedA> & A)
{
  assert(V.size() == S.size() && "S and V should be same size");
  A.setZero(S.maxCoeff()+1,1);
  for(int s = 0;s<S.size();s++)
  {
    A(S(s)) += V(s);
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiations
template void igl::accumarray<Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
#endif
