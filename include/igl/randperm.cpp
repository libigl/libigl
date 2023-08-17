// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "randperm.h"
#include "colon.h"
#include <algorithm>

template <typename DerivedI, typename URBG>
IGL_INLINE void igl::randperm(
  const int n,
  Eigen::PlainObjectBase<DerivedI> & I,
  URBG && urbg)
{
  Eigen::VectorXi II;
  igl::colon(0,1,n-1,II);
  I = II;

  std::shuffle(I.data(),I.data()+n, urbg);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::randperm<Eigen::Matrix<int, -1, -1, 0, -1, -1>, std::linear_congruential_engine<unsigned int, 48271u, 0u, 2147483647u>&>(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, std::linear_congruential_engine<unsigned int, 48271u, 0u, 2147483647u>&);
template void igl::randperm<Eigen::Matrix<int, -1, -1, 0, -1, -1>, std::linear_congruential_engine<unsigned int, 16807u, 0u, 2147483647u>&>(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, std::linear_congruential_engine<unsigned int, 16807u, 0u, 2147483647u>&);
template void igl::randperm<Eigen::Matrix<int, -1, -1, 0, -1, -1>, std::mersenne_twister_engine<unsigned long long, 64ul, 312ul, 156ul, 31ul, 13043109905998158313ull, 29ul, 6148914691236517205ull, 17ul, 8202884508482404352ull, 37ul, 18444473444759240704ull, 43ul, 6364136223846793005ull>&>(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, std::mersenne_twister_engine<unsigned long long, 64ul, 312ul, 156ul, 31ul, 13043109905998158313ull, 29ul, 6148914691236517205ull, 17ul, 8202884508482404352ull, 37ul, 18444473444759240704ull, 43ul, 6364136223846793005ull>&);
template void igl::randperm<Eigen::Matrix<int, -1, -1, 0, -1, -1>, std::mersenne_twister_engine<unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615u, 11ul, 4294967295u, 7ul, 2636928640u, 15ul, 4022730752u, 18ul, 1812433253u>&>(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, std::mersenne_twister_engine<unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615u, 11ul, 4294967295u, 7ul, 2636928640u, 15ul, 4022730752u, 18ul, 1812433253u>&);
template void igl::randperm<Eigen::Matrix<int, -1, -1, 0, -1, -1>, std::mersenne_twister_engine<unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615u, 11ul, 4294967295u, 7ul, 2636928640u, 15ul, 4022730752u, 18ul, 1812433253u> >(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, std::mersenne_twister_engine<unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615u, 11ul, 4294967295u, 7ul, 2636928640u, 15ul, 4022730752u, 18ul, 1812433253u>&&);
template void igl::randperm<Eigen::Matrix<int, -1, 1, 0, -1, 1>, std::linear_congruential_engine<unsigned int, 48271u, 0u, 2147483647u>&>(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, std::linear_congruential_engine<unsigned int, 48271u, 0u, 2147483647u>&);
template void igl::randperm<Eigen::Matrix<int, -1, 1, 0, -1, 1>, std::linear_congruential_engine<unsigned int, 16807u, 0u, 2147483647u>&>(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, std::linear_congruential_engine<unsigned int, 16807u, 0u, 2147483647u>&);
template void igl::randperm<Eigen::Matrix<int, -1, 1, 0, -1, 1>, std::mersenne_twister_engine<unsigned long long, 64ul, 312ul, 156ul, 31ul, 13043109905998158313ull, 29ul, 6148914691236517205ull, 17ul, 8202884508482404352ull, 37ul, 18444473444759240704ull, 43ul, 6364136223846793005ull>&>(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, std::mersenne_twister_engine<unsigned long long, 64ul, 312ul, 156ul, 31ul, 13043109905998158313ull, 29ul, 6148914691236517205ull, 17ul, 8202884508482404352ull, 37ul, 18444473444759240704ull, 43ul, 6364136223846793005ull>&);
template void igl::randperm<Eigen::Matrix<int, -1, 1, 0, -1, 1>, std::mersenne_twister_engine<unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615u, 11ul, 4294967295u, 7ul, 2636928640u, 15ul, 4022730752u, 18ul, 1812433253u>&>(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, std::mersenne_twister_engine<unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615u, 11ul, 4294967295u, 7ul, 2636928640u, 15ul, 4022730752u, 18ul, 1812433253u>&);
template void igl::randperm<Eigen::Matrix<int, -1, 1, 0, -1, 1>, std::mersenne_twister_engine<unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615u, 11ul, 4294967295u, 7ul, 2636928640u, 15ul, 4022730752u, 18ul, 1812433253u> >(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, std::mersenne_twister_engine<unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615u, 11ul, 4294967295u, 7ul, 2636928640u, 15ul, 4022730752u, 18ul, 1812433253u>&&);
// For some reason linux is using unsigned long 
template void igl::randperm<Eigen::Matrix<int, -1, -1, 0, -1, -1>, std::linear_congruential_engine<unsigned long, 16807ul, 0ul, 2147483647ul>&>(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, std::linear_congruential_engine<unsigned long, 16807ul, 0ul, 2147483647ul>&);
template void igl::randperm<Eigen::Matrix<int, -1, -1, 0, -1, -1>, std::linear_congruential_engine<unsigned long, 48271ul, 0ul, 2147483647ul>&>(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, std::linear_congruential_engine<unsigned long, 48271ul, 0ul, 2147483647ul>&);
template void igl::randperm<Eigen::Matrix<int, -1, -1, 0, -1, -1>, std::mersenne_twister_engine<unsigned long, 32ul, 624ul, 397ul, 31ul, 2567483615ul, 11ul, 4294967295ul, 7ul, 2636928640ul, 15ul, 4022730752ul, 18ul, 1812433253ul> >(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, std::mersenne_twister_engine<unsigned long, 32ul, 624ul, 397ul, 31ul, 2567483615ul, 11ul, 4294967295ul, 7ul, 2636928640ul, 15ul, 4022730752ul, 18ul, 1812433253ul>&&);
template void igl::randperm<Eigen::Matrix<int, -1, -1, 0, -1, -1>, std::mersenne_twister_engine<unsigned long, 32ul, 624ul, 397ul, 31ul, 2567483615ul, 11ul, 4294967295ul, 7ul, 2636928640ul, 15ul, 4022730752ul, 18ul, 1812433253ul>&>(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, std::mersenne_twister_engine<unsigned long, 32ul, 624ul, 397ul, 31ul, 2567483615ul, 11ul, 4294967295ul, 7ul, 2636928640ul, 15ul, 4022730752ul, 18ul, 1812433253ul>&);
template void igl::randperm<Eigen::Matrix<int, -1, -1, 0, -1, -1>, std::mersenne_twister_engine<unsigned long, 64ul, 312ul, 156ul, 31ul, 13043109905998158313ul, 29ul, 6148914691236517205ul, 17ul, 8202884508482404352ul, 37ul, 18444473444759240704ul, 43ul, 6364136223846793005ul>&>(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, std::mersenne_twister_engine<unsigned long, 64ul, 312ul, 156ul, 31ul, 13043109905998158313ul, 29ul, 6148914691236517205ul, 17ul, 8202884508482404352ul, 37ul, 18444473444759240704ul, 43ul, 6364136223846793005ul>&);
template void igl::randperm<Eigen::Matrix<int, -1, 1, 0, -1, 1>, std::linear_congruential_engine<unsigned long, 16807ul, 0ul, 2147483647ul>&>(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, std::linear_congruential_engine<unsigned long, 16807ul, 0ul, 2147483647ul>&);
template void igl::randperm<Eigen::Matrix<int, -1, 1, 0, -1, 1>, std::linear_congruential_engine<unsigned long, 48271ul, 0ul, 2147483647ul>&>(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, std::linear_congruential_engine<unsigned long, 48271ul, 0ul, 2147483647ul>&);
template void igl::randperm<Eigen::Matrix<int, -1, 1, 0, -1, 1>, std::mersenne_twister_engine<unsigned long, 32ul, 624ul, 397ul, 31ul, 2567483615ul, 11ul, 4294967295ul, 7ul, 2636928640ul, 15ul, 4022730752ul, 18ul, 1812433253ul> >(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, std::mersenne_twister_engine<unsigned long, 32ul, 624ul, 397ul, 31ul, 2567483615ul, 11ul, 4294967295ul, 7ul, 2636928640ul, 15ul, 4022730752ul, 18ul, 1812433253ul>&&);
template void igl::randperm<Eigen::Matrix<int, -1, 1, 0, -1, 1>, std::mersenne_twister_engine<unsigned long, 32ul, 624ul, 397ul, 31ul, 2567483615ul, 11ul, 4294967295ul, 7ul, 2636928640ul, 15ul, 4022730752ul, 18ul, 1812433253ul>&>(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, std::mersenne_twister_engine<unsigned long, 32ul, 624ul, 397ul, 31ul, 2567483615ul, 11ul, 4294967295ul, 7ul, 2636928640ul, 15ul, 4022730752ul, 18ul, 1812433253ul>&);
template void igl::randperm<Eigen::Matrix<int, -1, 1, 0, -1, 1>, std::mersenne_twister_engine<unsigned long, 64ul, 312ul, 156ul, 31ul, 13043109905998158313ul, 29ul, 6148914691236517205ul, 17ul, 8202884508482404352ul, 37ul, 18444473444759240704ul, 43ul, 6364136223846793005ul>&>(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, std::mersenne_twister_engine<unsigned long, 64ul, 312ul, 156ul, 31ul, 13043109905998158313ul, 29ul, 6148914691236517205ul, 17ul, 8202884508482404352ul, 37ul, 18444473444759240704ul, 43ul, 6364136223846793005ul>&);
#endif
