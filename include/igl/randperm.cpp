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

template <typename DerivedI>
IGL_INLINE void igl::randperm(
  const int n,
  Eigen::PlainObjectBase<DerivedI> & I,
  const int64_t rng_min,
  const int64_t rng_max,
  const std::function<int64_t()> &rng)
{
  Eigen::VectorXi II;
  igl::colon(0,1,n-1,II);
  I = II;

  // C++ named requirement : UniformRandomBitGenerator
  // This signature is required for the third parameter of
  // std::shuffle
  struct RandPermURBG {
  public:
    using result_type = int64_t;
    RandPermURBG(const result_type min,
                 const result_type max,
                 const std::function<result_type()> int_gen) :
      m_min(min),
      m_max(max),
      m_int_gen(std::move(int_gen)) {}

    result_type min() const noexcept { return m_min; }
    result_type max() const noexcept { return m_max; }
    result_type operator()() const { return m_int_gen(); };
  private:
    result_type m_min;
    result_type m_max;
    std::function<result_type()> m_int_gen;
  };

  const auto int_gen = (nullptr != rng) ?
    rng :
    []()->int64_t { return std::rand(); };

  std::shuffle(I.data(),I.data()+n, RandPermURBG(rng_min, rng_max, int_gen));
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::randperm<Eigen::Matrix<int, -1, 1, 0, -1, 1> >(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, const int64_t, const int64_t, const std::function<int64_t()>&);
template void igl::randperm<Eigen::Matrix<int, -1, -1, 0, -1, -1> >(int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, const int64_t, const int64_t, const std::function<int64_t()>&);
#endif
