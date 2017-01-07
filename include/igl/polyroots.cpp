// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "polyroots.h"
#include <Eigen/Eigenvalues>

template <typename S, typename T>
IGL_INLINE void igl::polyRoots(Eigen::Matrix<S, Eigen::Dynamic,1> &polyCoeff, //real or comples coefficients
                          Eigen::Matrix<std::complex<T>, Eigen::Dynamic,1> &roots // complex roots (double or float)
)
{
  //  degree
  int n = polyCoeff.rows() - 1;

  Eigen::Matrix<S, Eigen::Dynamic, 1> d (n,1);
  d = polyCoeff.tail(n)/polyCoeff(0);

  Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic> I; I.setIdentity(n-1,n-1);
  Eigen::Matrix<S, Eigen::Dynamic, 1> z; z.setZero(n-1,1);

  Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic> a(n,n);
  a<<-d.transpose(),I,z;
  roots = a.eigenvalues();

}



#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::polyRoots<std::complex<double>, double>(Eigen::Matrix<std::complex<double>, -1, 1, 0, -1, 1>&, Eigen::Matrix<std::complex<double>, -1, 1, 0, -1, 1>&);
template void igl::polyRoots<double, double>(Eigen::Matrix<double, -1, 1, 0, -1, 1>&, Eigen::Matrix<std::complex<double>, -1, 1, 0, -1, 1>&);
#endif
