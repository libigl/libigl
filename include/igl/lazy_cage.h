// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2024 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_LAZY_CAGE_H
#define IGL_LAZY_CAGE_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  /// Compute a "lazy cage" enclosing a closed input mesh: the coarse cage is
  /// obtained by offsetting (V,F) outward by the _smallest_ amount σ for which
  /// the dense offset surface can be decimated down to `num_faces` faces without
  /// ever intersecting the input surface or itself.
  ///
  /// Because collapses that would cross the input or create self-intersections
  /// are refused, a valid decimation to `num_faces` implies the cage stays on
  /// the outside of (and therefore encloses) the input and is free of
  /// self-intersections. The offset amount is found by bisection on σ: too small
  /// an offset leaves no room to decimate to the target without hitting the
  /// input; a large enough offset does.
  ///
  /// @param[in] V  #V by 3 list of input vertex positions (closed, manifold)
  /// @param[in] F  #F by 3 list of input triangle indices into V
  /// @param[in] num_faces  desired number of faces in the cage
  /// @param[in] grid_size  isosurfacing resolution: number of cells along the
  ///   largest bounding-box side used to sample the offset (larger = denser
  ///   offset, more faithful but slower)
  /// @param[in] max_sigma  upper bound on the searched offset distance
  /// @param[in] num_iters  number of bisection iterations on σ
  /// @param[out] CV  #CV by 3 list of cage vertex positions
  /// @param[out] CF  #CF by 3 list of cage triangle indices into CV
  /// @param[out] sigma  the offset distance used to build the returned cage
  /// @return true if a cage with `num_faces` faces was found within
  ///   `[0,max_sigma]`, false otherwise (in which case CV,CF hold the best
  ///   attempt, if any, and `sigma` its offset)
  ///
  /// \see decimate, block_self_intersections, block_intersections_with_input,
  ///   signed_distance, marching_cubes
  IGL_INLINE bool lazy_cage(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const int num_faces,
    const int grid_size,
    const double max_sigma,
    const int num_iters,
    Eigen::MatrixXd & CV,
    Eigen::MatrixXi & CF,
    double & sigma);
  /// \overload
  /// \brief Uses max_sigma = 0.1*bbox-diagonal and num_iters = 10.
  IGL_INLINE bool lazy_cage(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const int num_faces,
    const int grid_size,
    Eigen::MatrixXd & CV,
    Eigen::MatrixXi & CF,
    double & sigma);
}
#ifndef IGL_STATIC_LIBRARY
#  include "lazy_cage.cpp"
#endif
#endif
