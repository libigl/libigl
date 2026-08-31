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
#include "block_self_intersections.h"
#include "block_intersections_with_input.h"
#include "bounding_box_diagonal.h"
#include <limits>
#include <vector>

IGL_INLINE bool igl::lazy_cage(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const int num_faces,
  const int grid_size,
  const double max_sigma,
  const int num_iters,
  Eigen::MatrixXd & CV,
  Eigen::MatrixXi & CF,
  double & sigma)
{
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

  // Attempt to build a cage at a given offset: extract the offset isosurface and
  // decimate it to `num_faces` while forbidding self-intersections and
  // intersections with the input. Returns whether the target was reached.
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
    return igl::decimate(OV,OF,num_faces,decorators,U,G,J,I);
  };

  // The upper bound must succeed; otherwise there is no cage in range.
  double hi = max_sigma;
  Eigen::MatrixXd U_hi;
  Eigen::MatrixXi G_hi;
  const bool hi_ok = try_sigma(hi, U_hi, G_hi);
  CV = U_hi; CF = G_hi; sigma = hi;
  if(!hi_ok) { return false; }

  // Bisection for the smallest σ that still reaches the target. `lo` is kept as
  // a known (assumed) failure and is never evaluated at 0 (where the offset
  // coincides with the input).
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
  return lazy_cage(V,F,num_faces,grid_size,0.1*diag,10,CV,CF,sigma);
}
