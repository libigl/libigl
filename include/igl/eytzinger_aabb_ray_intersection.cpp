// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "eytzinger_aabb_ray_intersection.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

template <
  bool FirstOnly,
  typename Derivedsource,
  typename Deriveddir,
  typename DerivedB,
  typename Derivedleaf>
IGL_INLINE void igl::eytzinger_aabb_ray_intersection(
  const Eigen::MatrixBase<Derivedsource> & source,
  const Eigen::MatrixBase<Deriveddir> & dir,
  const std::function<bool(const int)> & hit,
  const std::function<bool(const int, const int)> & before,
  const Eigen::MatrixBase<DerivedB> & B1,
  const Eigen::MatrixBase<DerivedB> & B2,
  const Eigen::MatrixBase<Derivedleaf> & leaf,
  std::vector<int> & hits)
{
  using Scalar = typename std::common_type<typename Derivedsource::Scalar,
    typename Deriveddir::Scalar,typename DerivedB::Scalar>::type;
  hits.clear();
  const int dim = source.size();
  if(dir.size()!=dim || B1.cols()!=dim || B2.cols()!=dim ||
     B1.rows()!=leaf.size() || B2.rows()!=leaf.size())
    throw std::invalid_argument("Incompatible ray tree dimensions");
  const Scalar inf = std::numeric_limits<Scalar>::infinity();
  // Spatial inflation and parametric slack make box culling conservative.
  const Scalar releps = Scalar(1024) * std::numeric_limits<Scalar>::epsilon();
  const Scalar tslack = Scalar(8) * std::numeric_limits<Scalar>::epsilon();
  const Scalar tiny = std::numeric_limits<Scalar>::min();

  // Conservative ray/box slab test.  Returns whether the ray may enter the
  // (slightly inflated) box [lo,hi] and, if so, a lower bound `tnear` on the ray
  // parameter of any point inside the box.  Inflation guarantees that a box the
  // ray genuinely grazes is never rejected.
  const auto ray_box = [&](const int i, Scalar & tnear) -> bool
  {
    Scalar t0 = Scalar(0);
    Scalar t1 = inf;
    for(int k = 0; k < dim; k++)
    {
      const Scalar lo = B1(i,k);
      const Scalar hi = B2(i,k);
      const Scalar o = (Scalar)source(k);
      const Scalar dk = (Scalar)dir(k);
      const Scalar ext = hi - lo;
      const Scalar eps =
        releps * (std::abs(lo) + std::abs(hi) + ext) + tiny;
      const Scalar L = lo - eps;
      const Scalar H = hi + eps;
      if(dk == Scalar(0))
      {
        if(o < L || o > H){ return false; }
      }else
      {
        const Scalar inv = Scalar(1) / dk;
        // A nonzero subnormal direction is not parallel. If its reciprocal or
        // a slab endpoint overflows, leave this axis unconstrained.
        if(!std::isfinite(inv)){ continue; }
        Scalar ta = (L - o) * inv;
        Scalar tb = (H - o) * inv;
        if(!std::isfinite(ta) || !std::isfinite(tb)){ continue; }
        if(ta > tb){ std::swap(ta, tb); }
        // Subtraction, reciprocal and multiplication each round once. Eight
        // epsilons cover their relative error and rounding this enlargement;
        // tiny covers absolute error from subnormal intermediates. Spatial
        // inflation alone cannot cover distant origins' parameter error.
        ta -= tslack * std::abs(ta) + tiny;
        tb += tslack * std::abs(tb) + tiny;
        if(ta > t0){ t0 = ta; }
        if(tb < t1){ t1 = tb; }
        if(t0 > t1){ return false; }
      }
    }
    tnear = t0;
    return true;
  };

  if(leaf.size() == 0){ return; }
  {
    Scalar tnear;
    if(!ray_box(0, tnear)){ return; }
  }

  int best = -1;

  // Depth-first traversal. Nearer child is visited first (a floating point
  // heuristic that never affects the result). Only the conservative slab test
  // skips subtrees; the current best hit does not prune the traversal.
  std::vector<int> stack;
  stack.reserve(64);
  stack.push_back(0);
  while(!stack.empty())
  {
    const int i = stack.back();
    stack.pop_back();

    const int l = (int)leaf(i);
    if(l >= 0)
    {
      if(hit(l))
      {
        if(FirstOnly)
        {
          if(best < 0 || before(l, best)){ best = l; }
        }else
        {
          hits.push_back(l);
        }
      }
      continue;
    }
    if(l == -2){ continue; }

    if(l != -1 || i > (leaf.size()-3)/2)
      throw std::invalid_argument("Incomplete ray tree internal node");

    const int left_i = 2*i + 1;
    const int right_i = 2*i + 2;
    if(right_i>=leaf.size())
      throw std::invalid_argument("Incomplete ray tree internal node");
    Scalar tl, tr;
    const bool ok_l = leaf(left_i)!=-2 && ray_box(left_i, tl);
    const bool ok_r = leaf(right_i)!=-2 && ray_box(right_i, tr);
    // Push so that the nearer child is processed first (popped last -> pushed last).
    if(ok_l && ok_r)
    {
      if(tl <= tr)
      {
        stack.push_back(right_i);
        stack.push_back(left_i);
      }else
      {
        stack.push_back(left_i);
        stack.push_back(right_i);
      }
    }else if(ok_l)
    {
      stack.push_back(left_i);
    }else if(ok_r)
    {
      stack.push_back(right_i);
    }
  }

  if(FirstOnly)
  {
    if(best >= 0){ hits.push_back(best); }
  }else
  {
    std::sort(hits.begin(), hits.end(),
      [&](const int a, const int b){ return before(a, b); });
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
#define IGL_EARI(FO,B) \
  template void igl::eytzinger_aabb_ray_intersection<FO, \
    Eigen::Matrix<double,1,3,1,1,3>, Eigen::Matrix<double,1,3,1,1,3>, \
    B, Eigen::Matrix<int,-1,1,0,-1,1>>( \
    const Eigen::MatrixBase<Eigen::Matrix<double,1,3,1,1,3>>&, \
    const Eigen::MatrixBase<Eigen::Matrix<double,1,3,1,1,3>>&, \
    const std::function<bool(const int)>&, \
    const std::function<bool(const int, const int)>&, \
    const Eigen::MatrixBase<B>&, \
    const Eigen::MatrixBase<B>&, \
    const Eigen::MatrixBase<Eigen::Matrix<int,-1,1,0,-1,1>>&, \
    std::vector<int>&);
using EariBox3 = Eigen::Matrix<double,-1,3>;
IGL_EARI(true,EariBox3)
IGL_EARI(false,EariBox3)
IGL_EARI(true,Eigen::MatrixXd)
IGL_EARI(false,Eigen::MatrixXd)
#undef IGL_EARI
#endif
