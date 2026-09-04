// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "exact_ray_mesh_intersect.h"
#include "exact_ray_triangle_intersect.h"
#include "exact_ray_triangle_hit_compare.h"
#include "eytzinger_aabb.h"
#include "eytzinger_aabb_ray_intersection.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <functional>

namespace igl { namespace exact_ray_mesh_intersect_detail
{
  // Ordinary (inexact) Möller–Trumbore solve for the intersection parameters
  // (t,u,v). Used only to build the approximate convenience Hit record and the
  // floating point parameter that drives the traversal's culling/ordering
  // heuristics — never the exact classification. The only guard is the
  // mathematically required division by an exact zero determinant (no tolerance):
  // on inputs the exact predicate has already classified as a clean hit this
  // essentially never triggers.
  inline bool approx_tuv(
    const Eigen::RowVector3d & O, const Eigen::RowVector3d & d,
    const Eigen::RowVector3d & A, const Eigen::RowVector3d & B,
    const Eigen::RowVector3d & C,
    double & t, double & u, double & v)
  {
    const Eigen::RowVector3d e1 = B - A;
    const Eigen::RowVector3d e2 = C - A;
    const Eigen::RowVector3d pv = d.cross(e2);
    const double det = e1.dot(pv);
    if(det == 0.0){ t = u = v = 0.0; return false; }
    const double inv = 1.0 / det;
    const Eigen::RowVector3d tv = O - A;
    u = tv.dot(pv) * inv;
    const Eigen::RowVector3d qv = tv.cross(e1);
    v = d.dot(qv) * inv;
    t = e2.dot(qv) * inv;
    return true;
  }

  // Build the per-face Eytzinger AABB tree and run the traversal, returning the
  // hit face indices (size <= 1 when FirstOnly).
  template <
    bool FirstOnly,
    typename Derivedsource,
    typename Deriveddir,
    typename DerivedV,
    typename DerivedF>
  inline void query(
    const Eigen::MatrixBase<Derivedsource> & source,
    const Eigen::MatrixBase<Deriveddir> & dir,
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    std::vector<int> & fids)
  {
    fids.clear();
    const int m = F.rows();
    if(m == 0){ return; }

    const Eigen::RowVector3d S((double)source(0),(double)source(1),(double)source(2));
    const Eigen::RowVector3d D((double)dir(0),(double)dir(1),(double)dir(2));

    // Per-face axis-aligned bounding boxes (in double, matching the exact leaf).
    Eigen::Matrix<double,Eigen::Dynamic,3> PB1(m,3), PB2(m,3);
    for(int f = 0; f < m; f++)
    {
      for(int c = 0; c < 3; c++)
      {
        double mn = (double)V(F(f,0),c);
        double mx = mn;
        for(int j = 1; j < 3; j++)
        {
          const double x = (double)V(F(f,j),c);
          mn = std::min(mn,x);
          mx = std::max(mx,x);
        }
        PB1(f,c) = mn; PB2(f,c) = mx;
      }
    }

    Eigen::Matrix<double,Eigen::Dynamic,3> B1, B2;
    Eigen::VectorXi leaf;
    igl::eytzinger_aabb(PB1, PB2, B1, B2, leaf);

    // Materialize fixed-size row vectors so the (statically instantiated) leaf
    // predicates and traversal are always called with concrete types rather than
    // Eigen expression templates.
    const auto row = [&V](const int i) -> Eigen::RowVector3d
    {
      return Eigen::RowVector3d((double)V(i,0),(double)V(i,1),(double)V(i,2));
    };
    const std::function<bool(const int)> hit =
      [&](const int f) -> bool
    {
      // Only clean transverse hits count. Coplanar configurations cannot be
      // ordered with these exact predicates (their intersection is a segment,
      // not a point), so they are treated as non-hits.
      return igl::exact_ray_triangle_intersect(
        S, D, row(F(f,0)), row(F(f,1)), row(F(f,2))) ==
        igl::RayTriangleIntersection::Hit;
    };
    const std::function<bool(const int, const int)> before =
      [&](const int f1, const int f2) -> bool
    {
      return igl::exact_ray_triangle_hit_compare(
        S, D,
        row(F(f1,0)), row(F(f1,1)), row(F(f1,2)),
        row(F(f2,0)), row(F(f2,1)), row(F(f2,2))) < 0;
    };

    // The closest hit and the sort order are decided by the exact `before`
    // comparator. First-hit uses no distance pruning: benchmarks showed the exact
    // box-entry-vs-best pruning predicate costs as much as it saves (the
    // conservative box cull already limits traversal to the ray's neighborhood),
    // so the simpler, pruning-free strategy is used.
    igl::eytzinger_aabb_ray_intersection<FirstOnly>(
      S, D, hit, before, B1, B2, leaf, fids);
  }

  // Construct an approximate hit record for face fid.
  template <typename Derivedsource, typename Deriveddir, typename DerivedV, typename DerivedF>
  inline igl::Hit<typename DerivedV::Scalar> approximate_hit(
    const Eigen::MatrixBase<Derivedsource> & source,
    const Eigen::MatrixBase<Deriveddir> & dir,
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    const int fid)
  {
    using Scalar = typename DerivedV::Scalar;
    const Eigen::RowVector3d A((double)V(F(fid,0),0),(double)V(F(fid,0),1),(double)V(F(fid,0),2));
    const Eigen::RowVector3d B((double)V(F(fid,1),0),(double)V(F(fid,1),1),(double)V(F(fid,1),2));
    const Eigen::RowVector3d C((double)V(F(fid,2),0),(double)V(F(fid,2),1),(double)V(F(fid,2),2));
    const Eigen::RowVector3d S((double)source(0),(double)source(1),(double)source(2));
    const Eigen::RowVector3d D((double)dir(0),(double)dir(1),(double)dir(2));
    double t = 0, u = 0, v = 0;
    approx_tuv(S, D, A, B, C, t, u, v);
    return igl::Hit<Scalar>{ fid, -1, (Scalar)u, (Scalar)v, (Scalar)t };
  }
}}

template <
  typename Derivedsource,
  typename Deriveddir,
  typename DerivedV,
  typename DerivedF>
IGL_INLINE bool igl::exact_ray_mesh_intersect(
  const Eigen::MatrixBase<Derivedsource> & source,
  const Eigen::MatrixBase<Deriveddir> & dir,
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  int & fid)
{
  std::vector<int> fids;
  igl::exact_ray_mesh_intersect_detail::query<true>(source, dir, V, F, fids);
  if(fids.empty()){ return false; }
  fid = fids.front();
  return true;
}

template <
  typename Derivedsource,
  typename Deriveddir,
  typename DerivedV,
  typename DerivedF>
IGL_INLINE bool igl::exact_ray_mesh_intersect_all(
  const Eigen::MatrixBase<Derivedsource> & source,
  const Eigen::MatrixBase<Deriveddir> & dir,
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  std::vector<int> & fids)
{
  igl::exact_ray_mesh_intersect_detail::query<false>(source, dir, V, F, fids);
  return !fids.empty();
}

template <
  typename Derivedsource,
  typename Deriveddir,
  typename DerivedV,
  typename DerivedF>
IGL_INLINE bool igl::exact_ray_mesh_intersect(
  const Eigen::MatrixBase<Derivedsource> & source,
  const Eigen::MatrixBase<Deriveddir> & dir,
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  igl::Hit<typename DerivedV::Scalar> & hit)
{
  int fid;
  if(!igl::exact_ray_mesh_intersect(source, dir, V, F, fid)){ return false; }
  hit = igl::exact_ray_mesh_intersect_detail::approximate_hit(source, dir, V, F, fid);
  return true;
}

template <
  typename Derivedsource,
  typename Deriveddir,
  typename DerivedV,
  typename DerivedF>
IGL_INLINE bool igl::exact_ray_mesh_intersect_all(
  const Eigen::MatrixBase<Derivedsource> & source,
  const Eigen::MatrixBase<Deriveddir> & dir,
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  std::vector<igl::Hit<typename DerivedV::Scalar>> & hits)
{
  std::vector<int> fids;
  igl::exact_ray_mesh_intersect_all(source, dir, V, F, fids);
  hits.clear();
  hits.reserve(fids.size());
  for(const int fid : fids)
  {
    hits.push_back(
      igl::exact_ray_mesh_intersect_detail::approximate_hit(source, dir, V, F, fid));
  }
  return !hits.empty();
}

#ifdef IGL_STATIC_LIBRARY
namespace { namespace ermi_types {
  using V3 = Eigen::Matrix<double,3,1,0,3,1>;    // Vector3d
  using R3 = Eigen::Matrix<double,1,3,1,1,3>;    // RowVector3d
  using MXd = Eigen::Matrix<double,-1,-1,0,-1,-1>;
  using MXi = Eigen::Matrix<int,-1,-1,0,-1,-1>;
}}
#define IGL_ERMI(S,D,V,F) \
  template bool igl::exact_ray_mesh_intersect<S,D,V,F>( \
    const Eigen::MatrixBase<S>&, const Eigen::MatrixBase<D>&, \
    const Eigen::MatrixBase<V>&, const Eigen::MatrixBase<F>&, int&); \
  template bool igl::exact_ray_mesh_intersect_all<S,D,V,F>( \
    const Eigen::MatrixBase<S>&, const Eigen::MatrixBase<D>&, \
    const Eigen::MatrixBase<V>&, const Eigen::MatrixBase<F>&, std::vector<int>&); \
  template bool igl::exact_ray_mesh_intersect<S,D,V,F>( \
    const Eigen::MatrixBase<S>&, const Eigen::MatrixBase<D>&, \
    const Eigen::MatrixBase<V>&, const Eigen::MatrixBase<F>&, igl::Hit<V::Scalar>&); \
  template bool igl::exact_ray_mesh_intersect_all<S,D,V,F>( \
    const Eigen::MatrixBase<S>&, const Eigen::MatrixBase<D>&, \
    const Eigen::MatrixBase<V>&, const Eigen::MatrixBase<F>&, std::vector<igl::Hit<V::Scalar>>&);
IGL_ERMI(ermi_types::V3, ermi_types::V3, ermi_types::MXd, ermi_types::MXi)
IGL_ERMI(ermi_types::R3, ermi_types::R3, ermi_types::MXd, ermi_types::MXi)
#undef IGL_ERMI
#endif
