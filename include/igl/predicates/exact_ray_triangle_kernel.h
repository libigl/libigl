// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PREDICATES_EXACT_RAY_TRIANGLE_KERNEL_H
#define IGL_PREDICATES_EXACT_RAY_TRIANGLE_KERNEL_H

#include <predicates.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace igl { namespace predicates { namespace exact_ray_triangle_kernel
{
  static_assert(std::is_same<IGL_PREDICATES_REAL, double>::value,
    "Exact ray predicates require the double-precision predicates library");

  inline int sign(const double value)
  {
    return (value > 0) - (value < 0);
  }

  // A column is either an original vector p or the exact difference p-q.
  struct Column
  {
    const double *p;
    const double *q;
  };

  inline int entry(const Column &column, const int axis, double *result)
  {
    if(!column.q)
    {
      result[0] = column.p[axis];
      return 1;
    }
    double value, tail;
    ::predicates_two_diff(column.p[axis], column.q[axis], &value, &tail);
    int length = 0;
    if(tail != 0) { result[length++] = tail; }
    if(value != 0 || length == 0) { result[length++] = value; }
    return length;
  }

  // A product of expansions of lengths m,n needs at most 2*m*n components.
  // Every buffer belongs to this invocation; no expansion escapes a predicate.
  template<int Capacity>
  inline int product(
    const int left_length, const double *left,
    const int right_length, const double *right, double *result)
  {
    double scratch[Capacity], scaled[Capacity];
    double *current = result, *next = scratch;
    int length = 1;
    current[0] = 0;
    for(int i = 0; i < right_length; ++i)
    {
      const int scaled_length = ::predicates_expansion_scale(
        left_length, left, right[i], scaled);
      if(i == 0)
      {
        std::copy_n(scaled, scaled_length, current);
        length = scaled_length;
      }else
      {
        length = ::predicates_expansion_sum(
          length, current, scaled_length, scaled, next);
        std::swap(current, next);
      }
    }
    if(current != result) { std::copy_n(current, length, result); }
    return length;
  }

  // Three differences have at most two components each. Each of the six
  // determinant monomials needs at most 2*(2*2*2)*2 = 32 components, hence 192
  // in total. A raw first column reduces the bound to 96.
  inline int det3_expansion(
    const Column &first, const Column &second, const Column &third,
    double *result)
  {
    const Column columns[3] = {first, second, third};
    double entries[3][3][2];
    int lengths[3][3];
    for(int axis = 0; axis < 3; ++axis)
    {
      for(int column = 0; column < 3; ++column)
      {
        lengths[axis][column] = entry(columns[column], axis, entries[axis][column]);
      }
    }
    const int permutations[6][3] = {{0,1,2}, {1,2,0}, {2,0,1}, {0,2,1}, {1,0,2}, {2,1,0}};
    double scratch[192], pair[8], term[32];
    double *current = result, *next = scratch;
    int length = 1;
    current[0] = 0;
    for(int permutation = 0; permutation < 6; ++permutation)
    {
      const int *axes = permutations[permutation];
      const int pair_length = product<8>(
        lengths[axes[0]][0], entries[axes[0]][0],
        lengths[axes[1]][1], entries[axes[1]][1], pair);
      const int term_length = product<32>(
        pair_length, pair, lengths[axes[2]][2], entries[axes[2]][2], term);
      if(permutation >= 3)
      {
        for(int i = 0; i < term_length; ++i) { term[i] = -term[i]; }
      }
      if(permutation == 0)
      {
        std::copy_n(term, term_length, current);
        length = term_length;
      }else
      {
        length = ::predicates_expansion_sum(length, current, term_length, term, next);
        std::swap(current, next);
      }
    }
    if(current != result) { std::copy_n(current, length, result); }
    return length;
  }

  // Return a certified nonzero sign, or zero when exact evaluation is needed.
  inline int det3_filter(const Column & c0, const Column & c1, const Column & c2)
  {
    double M[3][3];
    const Column * col[3] = {&c0, &c1, &c2};
    for(int c = 0; c < 3; c++)
    {
      for(int r = 0; r < 3; r++)
      {
        M[r][c] = col[c]->q ? (col[c]->p[r] - col[c]->q[r]) : col[c]->p[r];
      }
    }
    const double det =
        M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
      - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
      + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]);
    const double a[3][3] = {
      {std::fabs(M[0][0]), std::fabs(M[0][1]), std::fabs(M[0][2])},
      {std::fabs(M[1][0]), std::fabs(M[1][1]), std::fabs(M[1][2])},
      {std::fabs(M[2][0]), std::fabs(M[2][1]), std::fabs(M[2][2])}};
    const double perm =
        a[0][0]*(a[1][1]*a[2][2] + a[1][2]*a[2][1])
      + a[0][1]*(a[1][0]*a[2][2] + a[1][2]*a[2][0])
      + a[0][2]*(a[1][0]*a[2][1] + a[1][1]*a[2][0]);
    // Generously conservative dynamic filter. When it cannot certify the sign
    // we fall back to expansion evaluation, subject to the documented arithmetic
    // range assumptions.
    const double errbound = 256.0 * DBL_EPSILON * perm;
    if(det > errbound){ return 1; }
    if(det < -errbound){ return -1; }
    return 0;
  }

  inline int det3_sign(const Column & c0, const Column & c1, const Column & c2)
  {
    const int filtered = det3_filter(c0, c1, c2);
    if(filtered != 0) { return filtered; }
    double expansion[192];
    const int length = det3_expansion(c0, c1, c2, expansion);
    return sign(expansion[length - 1]);
  }

  inline int classify(const double *origin, const double *direction,
    const double *a, const double *b, const double *c)
  {
    const Column ray{direction, nullptr}, ab{b,a}, ac{c,a};
    const int denominator = det3_sign(ray, ab, ac);
    int numerator = det3_filter(ab, ac, {a,origin});
    if(numerator == 0)
    {
      // Shewchuk orient3d(a,b,c,origin) has the sign of -n.(origin-a).
      // Its magnitude is deliberately not used as an exact construction.
      double pa[3], pb[3], pc[3], po[3];
      std::copy_n(a, 3, pa); std::copy_n(b, 3, pb);
      std::copy_n(c, 3, pc); std::copy_n(origin, 3, po);
      numerator = sign(::orient3d(pa, pb, pc, po));
    }
    if(denominator == 0) { return numerator == 0 ? 2 : 0; }
    if(det3_sign(ray, {a,origin}, {b,origin}) == -denominator) { return 0; }
    if(det3_sign(ray, {b,origin}, {c,origin}) == -denominator) { return 0; }
    if(det3_sign(ray, {c,origin}, {a,origin}) == -denominator) { return 0; }
    return numerator == 0 || numerator == denominator ? 1 : 0;
  }

  inline void hit_fraction(const double *origin, const double *direction,
    const double *a, const double *b, const double *c,
    double *numerator, int &numerator_length, double *denominator, int &denominator_length)
  {
    numerator_length = det3_expansion({b,a}, {c,a}, {a,origin}, numerator);
    denominator_length = det3_expansion({direction,nullptr}, {b,a}, {c,a}, denominator);
    if(denominator[denominator_length - 1] < 0)
    {
      for(int i = 0; i < numerator_length; ++i) { numerator[i] = -numerator[i]; }
      for(int i = 0; i < denominator_length; ++i) { denominator[i] = -denominator[i]; }
    }
  }

  inline int compare(const double *origin, const double *direction,
    const double *a1, const double *b1, const double *c1,
    const double *a2, const double *b2, const double *c2)
  {
    double numerator1[192], denominator1[96], numerator2[192], denominator2[96];
    int numerator1_length, denominator1_length, numerator2_length, denominator2_length;
    hit_fraction(origin, direction, a1,b1,c1,
      numerator1, numerator1_length, denominator1, denominator1_length);
    hit_fraction(origin, direction, a2,b2,c2,
      numerator2, numerator2_length, denominator2, denominator2_length);

    // A finite nonoverlapping double expansion occupies disjoint bit positions
    // between 2^-1074 and 2^1023, so it has at most 2098 nonzero components.
    // A scaled numerator adds at most 384. Reserving their sum bounds every
    // write, without the 73728-component polynomial bound or heap allocations.
    constexpr int max_components = std::numeric_limits<double>::max_exponent
      - std::numeric_limits<double>::min_exponent + std::numeric_limits<double>::digits;
    double accumulator[max_components + 384], scratch[max_components + 384], scaled[384];
    double *current = accumulator, *next = scratch;
    int length = 1;
    current[0] = 0;
    for(int side = 0; side < 2; ++side)
    {
      const double *numerator = side == 0 ? numerator1 : numerator2;
      const double *denominator = side == 0 ? denominator2 : denominator1;
      const int numerator_length = side == 0 ? numerator1_length : numerator2_length;
      const int denominator_length = side == 0 ? denominator2_length : denominator1_length;
      for(int i = 0; i < denominator_length; ++i)
      {
        const int scaled_length = ::predicates_expansion_scale(numerator_length,
          numerator, side == 0 ? denominator[i] : -denominator[i], scaled);
        length = ::predicates_expansion_sum(length, current, scaled_length, scaled, next);
        if(length > max_components)
        {
          throw std::overflow_error("Exact ray comparison exceeded the finite expansion capacity");
        }
        std::swap(current, next);
      }
    }
    return sign(current[length - 1]);
  }
}}}
#endif
