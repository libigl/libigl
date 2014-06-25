// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "histc.h"
#include "matlab_format.h"
#include <cassert>
#include <iostream>

template <typename DerivedX, typename DerivedE, typename DerivedN, typename DerivedB>
IGL_INLINE void igl::histc(
  const Eigen::PlainObjectBase<DerivedX > & X,
  const Eigen::PlainObjectBase<DerivedE > & E,
  Eigen::PlainObjectBase<DerivedN > & N,
  Eigen::PlainObjectBase<DerivedB > & B)
{
  histc(X,E,B);
  const int n = E.size();
  const int m = X.size();
  assert(m == B.size());
  N.resize(n,1);
  N.setConstant(0);
#pragma omp parallel for
  for(int j = 0;j<m;j++)
  {
    if(B(j) >= 0)
    {
#pragma omp atomic
      N(B(j))++;
    }
  }
}

template <typename DerivedX, typename DerivedE, typename DerivedB>
IGL_INLINE void igl::histc(
  const Eigen::PlainObjectBase<DerivedX > & X,
  const Eigen::PlainObjectBase<DerivedE > & E,
  Eigen::PlainObjectBase<DerivedB > & B)
{
  const int m = X.size();
  using namespace std;
  assert( 
    (E.bottomRightCorner(E.size()-1,1) -
      E.topLeftCorner(E.size()-1,1)).maxCoeff() >= 0 && 
    "E should be monotonically increasing");
  B.resize(m,1);
#pragma omp parallel for
  for(int j = 0;j<m;j++)
  {
    const double x = X(j);
    // Boring one-offs
    if(x < E(0) || x > E(E.size()-1))
    {
      B(j) = -1;
      continue;
    }
    // Find x in E
    int l = 0;
    int h = E.size()-1;
    int k = l;
    while((h-l)>1)
    {
      assert(x >= E(l));
      assert(x <= E(h));
      k = (h+l)/2;
      if(x < E(k))
      {
        h = k;
      }else
      {
        l = k;
      }
    }
    if(x == E(h))
    {
      k = h;
    }else
    {
      k = l;
    }
    B(j) = k;
  }
}

template <typename DerivedE>
IGL_INLINE void igl::histc(
    const typename DerivedE::Scalar & x,
    const Eigen::PlainObjectBase<DerivedE > & E,
    typename DerivedE::Index & b)
{
  Eigen::Matrix<typename DerivedE::Scalar,1,1> X;
  X(0) = x;
  Eigen::Matrix<typename DerivedE::Index,1,1> B;
  hist(X,E,B);
  b = B(0);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::histc<Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
template void igl::histc<Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::histc<Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::histc<Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::histc<Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
