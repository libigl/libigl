// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2024 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "lazy_cage.h"
#include "signed_distance.h"
#include "voxel_grid.h"
#include "marching_cubes.h"
#include "decimate.h"
#include "qslim.h"
#include "block_self_intersections.h"
#include "block_intersections_with_input.h"
#include "bounding_box_diagonal.h"
#include "centroid.h"
#include "doublearea.h"
#include "AABB.h"
#include "tri_tri_intersect.h"
#include <Eigen/Geometry>
#include <cmath>
#include <limits>
#include <vector>

namespace
{
  // Does any face of (AV,AF) transversally (non-coplanar) intersect (BV,BF),
  // accelerated by a prebuilt AABB tree over (BV,BF)?
  bool meshes_intersect(
    const Eigen::MatrixXd & AV,
    const Eigen::MatrixXi & AF,
    const Eigen::MatrixXd & BV,
    const Eigen::MatrixXi & BF,
    const igl::AABB<Eigen::MatrixXd,3> & Btree)
  {
    for(int f = 0; f < AF.rows(); f++)
    {
      const Eigen::RowVector3d A0 = AV.row(AF(f,0));
      const Eigen::RowVector3d A1 = AV.row(AF(f,1));
      const Eigen::RowVector3d A2 = AV.row(AF(f,2));
      Eigen::AlignedBox<double,3> box;
      box.extend(A0.transpose()).extend(A1.transpose()).extend(A2.transpose());
      std::vector<const igl::AABB<Eigen::MatrixXd,3>*> cand;
      Btree.append_intersecting_leaves(box,cand);
      for(const auto * c : cand)
      {
        const int g = c->m_primitive;
        const Eigen::RowVector3d B0 = BV.row(BF(g,0));
        const Eigen::RowVector3d B1 = BV.row(BF(g,1));
        const Eigen::RowVector3d B2 = BV.row(BF(g,2));
        bool coplanar = false;
        Eigen::RowVector3d s,t;
        if(igl::tri_tri_intersection_test_3d(A0,A1,A2,B0,B1,B2,coplanar,s,t) &&
           !coplanar)
        {
          return true;
        }
      }
    }
    return false;
  }

  // Objective value of a valid cage (smaller is better) under `metric`.
  double cage_metric(
    const Eigen::MatrixXd & CV,
    const Eigen::MatrixXi & CF,
    const double s,
    const igl::LazyCageMetric metric)
  {
    switch(metric)
    {
      case igl::LAZY_CAGE_METRIC_VOLUME:
      {
        Eigen::RowVector3d c; double vol = 0;
        igl::centroid(CV,CF,c,vol);
        return std::abs(vol);
      }
      case igl::LAZY_CAGE_METRIC_SURFACE_AREA:
      {
        Eigen::VectorXd dblA;
        igl::doublearea(CV,CF,dblA);
        return 0.5*dblA.sum();
      }
      case igl::LAZY_CAGE_METRIC_SIGMA:
      default:
        return s;
    }
  }
}

IGL_INLINE bool igl::lazy_cage(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const int num_faces,
  const int grid_size,
  const double max_sigma,
  const int num_iters,
  const bool use_qslim,
  const LazyCageMetric metric,
  Eigen::MatrixXd & CV,
  Eigen::MatrixXi & CF,
  double & sigma)
{
  const double diag = igl::bounding_box_diagonal(V);

  // Sample a signed distance field on a background grid padded to comfortably
  // contain the largest offset we might try. The field is independent of σ, so
  // it is computed once and every candidate offset is just a different marching
  // cubes isolevel.
  Eigen::MatrixXd GV;
  Eigen::RowVector3i side;
  igl::voxel_grid(V, max_sigma*1.5, grid_size, 1, GV, side);
  Eigen::VectorXd S;
  {
    Eigen::VectorXi I; Eigen::MatrixXd C,N;
    igl::signed_distance(
      GV, V, F, igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER,
      -std::numeric_limits<double>::infinity(),
       std::numeric_limits<double>::infinity(), S, I, C, N);
  }

  // Tree over the input, reused to validate every candidate cage.
  igl::AABB<Eigen::MatrixXd,3> input_tree;
  input_tree.init(V,F);

  // Attempt to build a valid cage at a given offset: extract the offset
  // isosurface, decimate it to `num_faces` while forbidding self-intersections
  // and intersections with the input, then verify the result is actually a
  // valid cage (reached the target, does not intersect the input, and encloses
  // it). Merely reaching the target is not sufficient: if the offset is under
  // the grid resolution the isosurface itself may already cross the input, and
  // the intersection-blocking callbacks only prevent *new* intersections.
  const auto try_sigma = [&](
    const double s, Eigen::MatrixXd & U, Eigen::MatrixXi & G) -> bool
  {
    Eigen::MatrixXd OV;
    Eigen::MatrixXi OF;
    igl::marching_cubes(S, GV, side(0), side(1), side(2), s, OV, OF);
    if(OF.rows() <= num_faces) { return false; }
    std::vector<igl::decimate_pre_post_collapse_callbacks_decorator> decorators;
    decorators.push_back(igl::block_self_intersections());
    decorators.push_back(igl::block_intersections_with_input(V,F));
    Eigen::VectorXi J,I;
    const bool reached = use_qslim ?
      igl::qslim   (OV,OF,num_faces,decorators,U,G,J,I) :
      igl::decimate(OV,OF,num_faces,decorators,U,G,J,I);
    if(!reached || G.rows() == 0) { return false; }
    // The cage must not transversally intersect the input.
    if(meshes_intersect(U,G,V,F,input_tree)) { return false; }
    // ...and it must enclose the input: since the cage does not intersect the
    // input, the (connected) input is wholly inside or wholly outside; require
    // every input vertex to be strictly inside the cage.
    Eigen::VectorXd Sd; Eigen::VectorXi Id; Eigen::MatrixXd Cd,Nd;
    igl::signed_distance(
      V,U,G,igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER,
      -std::numeric_limits<double>::infinity(),
       std::numeric_limits<double>::infinity(),Sd,Id,Cd,Nd);
    if((Sd.array() > 1e-7*diag).count() > 0) { return false; }
    return true;
  };

  // The upper bound must succeed; otherwise there is no valid cage in range.
  double hi = max_sigma;
  Eigen::MatrixXd U_hi;
  Eigen::MatrixXi G_hi;
  const bool hi_ok = try_sigma(hi, U_hi, G_hi);
  CV = U_hi; CF = G_hi; sigma = hi;
  if(!hi_ok) { return false; }

  // Bisection for the smallest σ that still yields a valid cage. `lo` is kept as
  // a known (assumed) failure and is never evaluated at 0 (where the offset
  // coincides with the input). This is the answer for the σ metric and the
  // lower end of the sweep for the volume/area metrics.
  double lo = 0.0;
  for(int it = 0; it < num_iters; it++)
  {
    const double mid = 0.5*(lo+hi);
    Eigen::MatrixXd U;
    Eigen::MatrixXi G;
    if(try_sigma(mid, U, G))
    {
      hi = mid;
      CV = U; CF = G; sigma = mid;
    }else
    {
      lo = mid;
    }
  }
  if(metric == LAZY_CAGE_METRIC_SIGMA) { return true; }

  // Volume / surface area are not monotonic in σ: giving decimation more room
  // (a larger offset) can remove more volume/area than the tightest offset
  // allows. Sweep σ from the valid floor up to max_sigma and keep the best.
  double best = cage_metric(CV,CF,sigma,metric);
  const double floor_sigma = sigma;
  for(int i = 1; i <= num_iters; i++)
  {
    const double s = floor_sigma + (max_sigma-floor_sigma)*(double(i)/num_iters);
    Eigen::MatrixXd U;
    Eigen::MatrixXi G;
    if(try_sigma(s, U, G))
    {
      const double m = cage_metric(U,G,s,metric);
      if(m < best)
      {
        best = m;
        CV = U; CF = G; sigma = s;
      }
    }
  }
  return true;
}

IGL_INLINE bool igl::lazy_cage(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const int num_faces,
  const int grid_size,
  Eigen::MatrixXd & CV,
  Eigen::MatrixXi & CF,
  double & sigma)
{
  const double diag = igl::bounding_box_diagonal(V);
  return lazy_cage(
    V,F,num_faces,grid_size,0.1*diag,12,false,
    LAZY_CAGE_METRIC_SIGMA,CV,CF,sigma);
}
