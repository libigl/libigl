// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_COLLAPSE_EDGE_H
#define IGL_COLLAPSE_EDGE_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>
#include <set>
namespace igl
{
  // Assumes (V,F) is a closed manifold mesh (except for previouslly collapsed
  // faces which should be set to: 
  // [IGL_COLLAPSE_EDGE_NULL IGL_COLLAPSE_EDGE_NULL IGL_COLLAPSE_EDGE_NULL].
  // Collapses exactly two faces and exactly 3 edges from E (e and one side of
  // each face gets collapsed to the other). This is implemented in a way that
  // it can be repeatedly called until satisfaction and then the garbage in F
  // can be collected by removing NULL faces.
  //
  // Inputs:
  //   e  index into E of edge to try to collapse. E(e,:) = [s d] or [d s] so
  //     that s<d, then d is collapsed to s.
  ///  p  dim list of vertex position where to place merged vertex
  // Inputs/Outputs:
  //   V  #V by dim list of vertex positions, lesser index of E(e,:) will be set
  //     to midpoint of edge.
  //   F  #F by 3 list of face indices into V.
  //   E  #E by 2 list of edge indices into V.
  //   EMAP #F*3 list of indices into E, mapping each directed edge to unique
  //     unique edge in E
  //   EF  #E by 2 list of edge flaps, EF(e,0)=f means e=(i-->j) is the edge of
  //     F(f,:) opposite the vth corner, where EI(e,0)=v. Similarly EF(e,1) "
  //     e=(j->i)
  //   EI  #E by 2 list of edge flap corners (see above).
  //   e1  index into E of edge collpased on left
  //   e2  index into E of edge collpased on left
  //   f1  index into E of edge collpased on left
  //   f2  index into E of edge collpased on left
  // Returns true if edge was collapsed
  #define IGL_COLLAPSE_EDGE_NULL 0
  IGL_INLINE bool collapse_edge(
    const int e,
    const Eigen::RowVectorXd & p,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & E,
    Eigen::VectorXi & EMAP,
    Eigen::MatrixXi & EF,
    Eigen::MatrixXi & EI,
    int & e1,
    int & e2,
    int & f1,
    int & f2);
  IGL_INLINE bool collapse_edge(
    const int e,
    const Eigen::RowVectorXd & p,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & E,
    Eigen::VectorXi & EMAP,
    Eigen::MatrixXi & EF,
    Eigen::MatrixXi & EI);
  // Collapse least-cost edge from a priority queue and update queue 
  //
  // Inputs/Outputs:
  //   cost_and_placement  function computing cost of collapsing an edge and 3d
  //     position where it should be placed:
  //     cost_and_placement(V,F,E,EMAP,EF,EI,cost,placement);
  //   Q  queue containing pairs of costs and edge indices
  //   Qit  list of iterators so that Qit[e] --> iterator of edge e in Q
  //   C  #E by dim list of stored placements
  IGL_INLINE bool collapse_edge(
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
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & E,
    Eigen::VectorXi & EMAP,
    Eigen::MatrixXi & EF,
    Eigen::MatrixXi & EI,
    std::set<std::pair<double,int> > & Q,
    std::vector<std::set<std::pair<double,int> >::iterator > & Qit,
    Eigen::MatrixXd & C);
  IGL_INLINE bool collapse_edge(
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
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & E,
    Eigen::VectorXi & EMAP,
    Eigen::MatrixXi & EF,
    Eigen::MatrixXi & EI,
    std::set<std::pair<double,int> > & Q,
    std::vector<std::set<std::pair<double,int> >::iterator > & Qit,
    Eigen::MatrixXd & C,
    int & e,
    int & e1,
    int & e2,
    int & f1,
    int & f2);
}

#ifndef IGL_STATIC_LIBRARY
#  include "collapse_edge.cpp"
#endif
#endif
