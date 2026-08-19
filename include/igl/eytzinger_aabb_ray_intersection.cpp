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

template <
  bool FirstOnly,
  typename Derivedsource,
  typename Deriveddir,
  typename DerivedB,
  typename Derivedleaf>
IGL_INLINE void igl::eytzinger_aabb_ray_intersection(
  const Eigen::MatrixBase<Derivedsource> & source,
  const Eigen::MatrixBase<Deriveddir> & dir,
  const std::function<bool(const int, double &)> & hit,
  const std::function<bool(const int, const int)> & before,
  const Eigen::MatrixBase<DerivedB> & B1,
  const Eigen::MatrixBase<DerivedB> & B2,
  const Eigen::MatrixBase<Derivedleaf> & leaf,
  std::vector<int> & hits)
{
  hits.clear();
  const int dim = source.size();
  const double inf = std::numeric_limits<double>::infinity();

  // Conservative ray/box slab test.  Returns whether the ray may enter the
  // (slightly inflated) box [lo,hi] and, if so, a lower bound `tnear` on the ray
  // parameter of any point inside the box.  Inflation guarantees that a box the
  // ray genuinely grazes is never rejected.
  const auto ray_box = [&](const int i, double & tnear) -> bool
  {
    double t0 = 0.0;
    double t1 = inf;
    for(int k = 0; k < dim; k++)
    {
      const double lo = (double)B1(i,k);
      const double hi = (double)B2(i,k);
      const double o = (double)source(k);
      const double dk = (double)dir(k);
      const double ext = hi - lo;
      const double eps =
        1e-12 * (std::fabs(lo) + std::fabs(hi) + ext) + 1e-300;
      const double L = lo - eps;
      const double H = hi + eps;
      if(std::fabs(dk) <= 1e-300)
      {
        if(o < L || o > H){ return false; }
      }else
      {
        const double inv = 1.0 / dk;
        double ta = (L - o) * inv;
        double tb = (H - o) * inv;
        if(ta > tb){ std::swap(ta, tb); }
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
    double tnear;
    if(!ray_box(0, tnear)){ return; }
  }

  int best = -1;
  double best_t = inf;

  // Depth-first traversal.  Entries carry the box's conservative entry distance so
  // subtrees that certainly start beyond the current best can be pruned.
  std::vector<std::pair<int,double>> stack;
  stack.reserve(64);
  stack.emplace_back(0, 0.0);
  while(!stack.empty())
  {
    const auto [i, tnear] = stack.back();
    stack.pop_back();

    if(FirstOnly && best >= 0)
    {
      // Safety margin absorbs the floating point error in best_t (an approximate
      // parameter of an actual hit) and in tnear.
      const double margin = 1e-6 * (std::fabs(best_t) + 1.0);
      if(tnear > best_t + margin){ continue; }
    }

    const int l = (int)leaf(i);
    if(l >= 0)
    {
      double t;
      if(hit(l, t))
      {
        if(FirstOnly)
        {
          if(best < 0 || before(l, best))
          {
            best = l;
            best_t = t;
          }
        }else
        {
          hits.push_back(l);
        }
      }
      continue;
    }
    if(l == -2){ continue; }

    const int left_i = 2*i + 1;
    const int right_i = 2*i + 2;
    double tl, tr;
    const bool ok_l = (left_i  < leaf.size()) && (leaf(left_i)  != -2) && ray_box(left_i,  tl);
    const bool ok_r = (right_i < leaf.size()) && (leaf(right_i) != -2) && ray_box(right_i, tr);
    // Push so that the nearer child is processed first (popped last -> pushed last).
    if(ok_l && ok_r)
    {
      if(tl <= tr)
      {
        stack.emplace_back(right_i, tr);
        stack.emplace_back(left_i,  tl);
      }else
      {
        stack.emplace_back(left_i,  tl);
        stack.emplace_back(right_i, tr);
      }
    }else if(ok_l)
    {
      stack.emplace_back(left_i, tl);
    }else if(ok_r)
    {
      stack.emplace_back(right_i, tr);
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
#define IGL_EARI(FO) \
  template void igl::eytzinger_aabb_ray_intersection<FO, \
    Eigen::Matrix<double,1,3,1,1,3>, Eigen::Matrix<double,1,3,1,1,3>, \
    Eigen::Matrix<double,-1,3,0,-1,3>, Eigen::Matrix<int,-1,1,0,-1,1>>( \
    const Eigen::MatrixBase<Eigen::Matrix<double,1,3,1,1,3>>&, \
    const Eigen::MatrixBase<Eigen::Matrix<double,1,3,1,1,3>>&, \
    const std::function<bool(const int, double&)>&, \
    const std::function<bool(const int, const int)>&, \
    const Eigen::MatrixBase<Eigen::Matrix<double,-1,3,0,-1,3>>&, \
    const Eigen::MatrixBase<Eigen::Matrix<double,-1,3,0,-1,3>>&, \
    const Eigen::MatrixBase<Eigen::Matrix<int,-1,1,0,-1,1>>&, \
    std::vector<int>&);
IGL_EARI(true)
IGL_EARI(false)
#undef IGL_EARI
#endif
