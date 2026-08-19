// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "ray_mesh_intersect.h"
#include "ray_triangle_intersect.h"
#include "../eytzinger_aabb.h"
#include "../eytzinger_aabb_ray_intersection.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <limits>

namespace igl { namespace predicates { namespace ray_mesh_intersect_detail
{
  // Approximate (inexact) Moller-Trumbore. Only used for culling/ordering and for
  // the convenience hit records. Never affects the exact classification.
  template <typename DerivedS, typename DerivedD, typename DerivedA>
  inline bool approx_hit(
    const Eigen::MatrixBase<DerivedS> & source,
    const Eigen::MatrixBase<DerivedD> & dir,
    const Eigen::MatrixBase<DerivedA> & A,
    const Eigen::MatrixBase<DerivedA> & B,
    const Eigen::MatrixBase<DerivedA> & C,
    double & t, double & u, double & v)
  {
    const Eigen::Vector3d O(source(0), source(1), source(2));
    const Eigen::Vector3d d(dir(0), dir(1), dir(2));
    const Eigen::Vector3d a(A(0), A(1), A(2));
    const Eigen::Vector3d b(B(0), B(1), B(2));
    const Eigen::Vector3d c(C(0), C(1), C(2));
    const Eigen::Vector3d e1 = b - a;
    const Eigen::Vector3d e2 = c - a;
    const Eigen::Vector3d pv = d.cross(e2);
    const double det = e1.dot(pv);
    if(det == 0.0){ t = 0.0; u = 0.0; v = 0.0; return false; }
    const double inv = 1.0 / det;
    const Eigen::Vector3d tv = O - a;
    u = tv.dot(pv) * inv;
    const Eigen::Vector3d qv = tv.cross(e1);
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

    // Per-face axis-aligned bounding boxes.
    Eigen::Matrix<double,Eigen::Dynamic,3> PB1(m,3), PB2(m,3);
    for(int f = 0; f < m; f++)
    {
      Eigen::RowVector3d lo, hi;
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
        lo(c) = mn; hi(c) = mx;
      }
      PB1.row(f) = lo; PB2.row(f) = hi;
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
    const std::function<bool(const int, double &)> hit =
      [&](const int f, double & t) -> bool
    {
      const Eigen::RowVector3d A = row(F(f,0));
      const Eigen::RowVector3d B = row(F(f,1));
      const Eigen::RowVector3d C = row(F(f,2));
      if(!igl::predicates::ray_triangle_intersect(S, D, A, B, C))
      {
        return false;
      }
      double u, v;
      if(!approx_hit(S, D, A, B, C, t, u, v))
      {
        // Exact predicate certifies a transverse hit but the (inexact) parameter
        // could not be formed (near-parallel in double). Use +inf so this leaf
        // never tightens the pruning bound; exact `before` still orders it.
        t = std::numeric_limits<double>::infinity();
      }
      return true;
    };
    const std::function<bool(const int, const int)> before =
      [&](const int f1, const int f2) -> bool
    {
      return igl::predicates::ray_triangle_hit_compare(
        S, D,
        row(F(f1,0)), row(F(f1,1)), row(F(f1,2)),
        row(F(f2,0)), row(F(f2,1)), row(F(f2,2))) < 0;
    };

    igl::eytzinger_aabb_ray_intersection<FirstOnly>(
      S, D, hit, before, B1, B2, leaf, fids);
  }
}}}

template <
  typename Derivedsource,
  typename Deriveddir,
  typename DerivedV,
  typename DerivedF>
IGL_INLINE bool igl::predicates::ray_mesh_intersect(
  const Eigen::MatrixBase<Derivedsource> & source,
  const Eigen::MatrixBase<Deriveddir> & dir,
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  int & fid)
{
  std::vector<int> fids;
  igl::predicates::ray_mesh_intersect_detail::query<true>(source, dir, V, F, fids);
  if(fids.empty()){ return false; }
  fid = fids.front();
  return true;
}

template <
  typename Derivedsource,
  typename Deriveddir,
  typename DerivedV,
  typename DerivedF>
IGL_INLINE bool igl::predicates::ray_mesh_intersect_all(
  const Eigen::MatrixBase<Derivedsource> & source,
  const Eigen::MatrixBase<Deriveddir> & dir,
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  std::vector<int> & fids)
{
  igl::predicates::ray_mesh_intersect_detail::query<false>(source, dir, V, F, fids);
  return !fids.empty();
}

template <
  typename Derivedsource,
  typename Deriveddir,
  typename DerivedV,
  typename DerivedF>
IGL_INLINE bool igl::predicates::ray_mesh_intersect(
  const Eigen::MatrixBase<Derivedsource> & source,
  const Eigen::MatrixBase<Deriveddir> & dir,
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  igl::Hit<typename DerivedV::Scalar> & hit)
{
  using Scalar = typename DerivedV::Scalar;
  int fid;
  if(!igl::predicates::ray_mesh_intersect(source, dir, V, F, fid)){ return false; }
  double t, u, v;
  igl::predicates::ray_mesh_intersect_detail::approx_hit(
    source, dir, V.row(F(fid,0)), V.row(F(fid,1)), V.row(F(fid,2)), t, u, v);
  hit = { fid, -1, (Scalar)u, (Scalar)v, (Scalar)t };
  return true;
}

template <
  typename Derivedsource,
  typename Deriveddir,
  typename DerivedV,
  typename DerivedF>
IGL_INLINE bool igl::predicates::ray_mesh_intersect_all(
  const Eigen::MatrixBase<Derivedsource> & source,
  const Eigen::MatrixBase<Deriveddir> & dir,
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  std::vector<igl::Hit<typename DerivedV::Scalar>> & hits)
{
  using Scalar = typename DerivedV::Scalar;
  std::vector<int> fids;
  igl::predicates::ray_mesh_intersect_all(source, dir, V, F, fids);
  hits.clear();
  hits.reserve(fids.size());
  for(const int fid : fids)
  {
    double t, u, v;
    igl::predicates::ray_mesh_intersect_detail::approx_hit(
      source, dir, V.row(F(fid,0)), V.row(F(fid,1)), V.row(F(fid,2)), t, u, v);
    hits.push_back({ fid, -1, (Scalar)u, (Scalar)v, (Scalar)t });
  }
  return !hits.empty();
}

#ifdef IGL_STATIC_LIBRARY
namespace { namespace rmi_types {
  using V3 = Eigen::Matrix<double,3,1,0,3,1>;    // Vector3d
  using R3 = Eigen::Matrix<double,1,3,1,1,3>;    // RowVector3d
  using MXd = Eigen::Matrix<double,-1,-1,0,-1,-1>;
  using MXi = Eigen::Matrix<int,-1,-1,0,-1,-1>;
}}
#define IGL_RMI(S,D,V,F) \
  template bool igl::predicates::ray_mesh_intersect<S,D,V,F>( \
    const Eigen::MatrixBase<S>&, const Eigen::MatrixBase<D>&, \
    const Eigen::MatrixBase<V>&, const Eigen::MatrixBase<F>&, int&); \
  template bool igl::predicates::ray_mesh_intersect_all<S,D,V,F>( \
    const Eigen::MatrixBase<S>&, const Eigen::MatrixBase<D>&, \
    const Eigen::MatrixBase<V>&, const Eigen::MatrixBase<F>&, std::vector<int>&); \
  template bool igl::predicates::ray_mesh_intersect<S,D,V,F>( \
    const Eigen::MatrixBase<S>&, const Eigen::MatrixBase<D>&, \
    const Eigen::MatrixBase<V>&, const Eigen::MatrixBase<F>&, igl::Hit<V::Scalar>&); \
  template bool igl::predicates::ray_mesh_intersect_all<S,D,V,F>( \
    const Eigen::MatrixBase<S>&, const Eigen::MatrixBase<D>&, \
    const Eigen::MatrixBase<V>&, const Eigen::MatrixBase<F>&, std::vector<igl::Hit<V::Scalar>>&);
IGL_RMI(rmi_types::V3, rmi_types::V3, rmi_types::MXd, rmi_types::MXi)
IGL_RMI(rmi_types::R3, rmi_types::R3, rmi_types::MXd, rmi_types::MXi)
#undef IGL_RMI
#endif
