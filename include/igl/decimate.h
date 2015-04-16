// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_DECIMATE_H
#define IGL_DECIMATE_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>
#include <set>
namespace igl
{
  // Assumes (V,F) is a closed manifold mesh 
  // Collapses edges until desired number of faces is achieved. This uses
  // default edge cost and merged vertex placement functions {edge length, edge
  // midpoint}.
  //
  // Inputs:
  //   V  #V by dim list of vertex positions
  //   F  #F by 3 list of face indices into V.
  //   max_m  desired number of output faces
  // Outputs:
  //   U  #U by dim list of output vertex posistions (can be same ref as V)
  //   G  #G by 3 list of output face indices into U (can be same ref as G)
  // Returns true if m was reached (otherwise #G > m)
  IGL_INLINE bool decimate(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const size_t max_m,
    Eigen::MatrixXd & U,
    Eigen::MatrixXi & G);
  // Inputs:
  //   cost_and_placement  function computing cost of collapsing an edge and 3d
  //     position where it should be placed:
  //     cost_and_placement(V,F,E,EMAP,EF,EI,cost,placement);
  //   stopping_condition  function returning whether to stop collapsing edges
  //     based on current state. Guaranteed to be called after _successfully_
  //     collapsing edge e removing edges (e,e1,e2) and faces (f1,f2):
  //     bool should_stop =
  //       stopping_condition(V,F,E,EMAP,EF,EI,Q,Qit,C,e,e1,e2,f1,f2);
  IGL_INLINE bool decimate(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const std::function<void(
      const int,
      const Eigen::MatrixXd &,
      const Eigen::MatrixXi &,
      const Eigen::MatrixXi &,
      const Eigen::VectorXi &,
      const Eigen::MatrixXi &,
      const Eigen::MatrixXi &,
      double &,
      Eigen::RowVectorXd &)> & cost_and_placement,
    const std::function<bool(
      const Eigen::MatrixXd &,
      const Eigen::MatrixXi &,
      const Eigen::MatrixXi &,
      const Eigen::VectorXi &,
      const Eigen::MatrixXi &,
      const Eigen::MatrixXi &,
      const std::set<std::pair<double,int> > &,
      const std::vector<std::set<std::pair<double,int> >::iterator > &,
      const Eigen::MatrixXd &,
      const int,
      const int,
      const int,
      const int,
      const int)> & stopping_condition,
    Eigen::MatrixXd & U,
    Eigen::MatrixXi & G);

}

#ifndef IGL_STATIC_LIBRARY
#  include "decimate.cpp"
#endif
#endif

