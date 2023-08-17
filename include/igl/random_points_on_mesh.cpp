// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "random_points_on_mesh.h"
#include "random_points_on_mesh_intrinsic.h"
#include "doublearea.h"


template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedB,
  typename DerivedFI,
  typename DerivedX,
  typename URBG>
IGL_INLINE void igl::random_points_on_mesh(
  const int n,
  const Eigen::MatrixBase<DerivedV > & V,
  const Eigen::MatrixBase<DerivedF > & F,
  Eigen::PlainObjectBase<DerivedB > & B,
  Eigen::PlainObjectBase<DerivedFI > & FI,
  Eigen::PlainObjectBase<DerivedX> & X,
  URBG && urbg)
{
  using namespace Eigen;
  using namespace std;
  typedef typename DerivedV::Scalar Scalar;
  typedef Matrix<Scalar,Dynamic,1> VectorXs;
  VectorXs dblA;
  doublearea(V,F,dblA);
  random_points_on_mesh_intrinsic(n,dblA,B,FI,urbg);
  X = DerivedX::Zero(B.rows(),V.cols());
  for(int x = 0;x<B.rows();x++)
  {
    for(int b = 0;b<B.cols();b++)
    {
      auto fi = FI(x);
      auto f = F(fi, b);
      auto bval = B(x, b);
      X.row(x) += bval*V.row(f);
    }
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::random_points_on_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, std::linear_congruential_engine<unsigned int, 16807u, 0u, 2147483647u>&>(int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, std::linear_congruential_engine<unsigned int, 16807u, 0u, 2147483647u>&);
template void igl::random_points_on_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, std::linear_congruential_engine<unsigned int, 48271u, 0u, 2147483647u>&>(int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, std::linear_congruential_engine<unsigned int, 48271u, 0u, 2147483647u>&);
template void igl::random_points_on_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, std::mersenne_twister_engine<unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615u, 11ul, 4294967295u, 7ul, 2636928640u, 15ul, 4022730752u, 18ul, 1812433253u> >(int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, std::mersenne_twister_engine<unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615u, 11ul, 4294967295u, 7ul, 2636928640u, 15ul, 4022730752u, 18ul, 1812433253u>&&);
template void igl::random_points_on_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, std::mersenne_twister_engine<unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615u, 11ul, 4294967295u, 7ul, 2636928640u, 15ul, 4022730752u, 18ul, 1812433253u>&>(int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, std::mersenne_twister_engine<unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615u, 11ul, 4294967295u, 7ul, 2636928640u, 15ul, 4022730752u, 18ul, 1812433253u>&);
template void igl::random_points_on_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, std::mersenne_twister_engine<unsigned long long, 64ul, 312ul, 156ul, 31ul, 13043109905998158313ull, 29ul, 6148914691236517205ull, 17ul, 8202884508482404352ull, 37ul, 18444473444759240704ull, 43ul, 6364136223846793005ull>&>(int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, std::mersenne_twister_engine<unsigned long long, 64ul, 312ul, 156ul, 31ul, 13043109905998158313ull, 29ul, 6148914691236517205ull, 17ul, 8202884508482404352ull, 37ul, 18444473444759240704ull, 43ul, 6364136223846793005ull>&);
template void igl::random_points_on_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, std::mersenne_twister_engine<unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615u, 11ul, 4294967295u, 7ul, 2636928640u, 15ul, 4022730752u, 18ul, 1812433253u> >(int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, std::mersenne_twister_engine<unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615u, 11ul, 4294967295u, 7ul, 2636928640u, 15ul, 4022730752u, 18ul, 1812433253u>&&);
template void igl::random_points_on_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 3, 1, -1, 3>, std::mersenne_twister_engine<unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615u, 11ul, 4294967295u, 7ul, 2636928640u, 15ul, 4022730752u, 18ul, 1812433253u>&>(int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&, std::mersenne_twister_engine<unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615u, 11ul, 4294967295u, 7ul, 2636928640u, 15ul, 4022730752u, 18ul, 1812433253u>&);
#endif
