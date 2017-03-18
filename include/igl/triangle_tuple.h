// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Zhongshi Jiang <zhongshi@cims.nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_TRIANGLE_TUPLE_H
#define IGL_TRIANGLE_TUPLE_H

#include <igl/igl_inline.h>
#include <Eigen/Core>

namespace igl
{
  // triangle cell tuple - Fake half-edge for fast and easy navigation
  // on trianglerahedral meshes with triangle_triangle_adjacency.

  // on triangle meshes with vertex_triangle_adjacency and
  // triangle_triangle adjacency
  //
  // Note: this is different to classical Half Edge data structure.
  //    Instead, it follows cell-tuple in [Brisson, 1989]
  //    "Representing geometric structures in d dimensions: topology and order."
  //    This class can achieve local navigation similar to half edge in OpenMesh
  //    But the logic behind each atom operation is different.
  //    So this should be more properly called TriangleTupleIterator.
  //
  // Each tuple contains information of (face, edge, (from)vertex) associated
  // with a directed half edge and encoded by (face, edge \in {0,1,2}, bool
  // along), stands for the index of triangle, local index of edge, and whether
  // the half edge is along the natural positive direction of the edge.
  // There are two type of atom operation: switch/get detailed below.
  //
  // Templates:
  //    DerivedF/FF/FFi Matrix Type for F/FF/FFi. Can be implicitly captured.
  // Inputs:
  //    (fi, ei, along) as detailed above.
  //    F #F by 3 list of triangular faces
  //    FF #F by 3 list of triangle-triangle adjacency.
  //    FFi #F by 3 list of FF inverse. Refer to
  //        "triangle_triangle_adjacency.h"
  // Usages:
  //    triangle_tuple_get_{vert/edge/face} returns the
  //      {index of from-vertex / *local* edge index / index of triangle}
  //      associated with the represented half edge.
  //
  //    triangle_tuple_switch_{vert/edge/face} switches to another half edge
  //      that differes in only one v/e/f information.
  //
  //    triangle_tuple_next_in_one_ring switch to next half edge
  //      in one-ring of the from vertex: (switch_face + switch_edge)
  //      with proper handling of boundary.
  //
  //    triangle_tuple_is_on_boundary determines whether the current encoded
  //      half edge is on boundary.
  //
  //    triangle_tuples_equal is only a wrapper comparing the (f,e,a) tuple
  
  // Basic Switch: only change information for one face/edge/vert resp.
  template<typename DerivedF, typename DerivedFF,
           typename DerivedFFi>
  IGL_INLINE void triangle_tuple_switch_vert(
      int &fi, int &ei, bool &along,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedFF> &FF,
      const Eigen::MatrixBase<DerivedFFi> &FFi);

  template<typename DerivedF, typename DerivedFF,
           typename DerivedFFi>
  IGL_INLINE void triangle_tuple_switch_edge(
      int &fi, int &ei, bool &along,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedFF> &FF,
      const Eigen::MatrixBase<DerivedFFi> &FFi);

  template<typename DerivedF, typename DerivedFF,
           typename DerivedFFi> 
  IGL_INLINE void triangle_tuple_switch_face(
      int &fi, int &ei, bool &along,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedFF> &FF,
      const Eigen::MatrixBase<DerivedFFi> &FFi);

  //getters
  template<typename DerivedF, typename DerivedFF,
           typename DerivedFFi>
  IGL_INLINE int triangle_tuple_get_vert(
      const int &fi, const int &ei, const bool &along,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedFF> &FF,
      const Eigen::MatrixBase<DerivedFFi> &FFi);

  template<typename DerivedF,
           typename DerivedFF,
           typename DerivedFFi>
  IGL_INLINE int triangle_tuple_get_edge(
      const int &fi, const int &ei, const bool &along,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedFF> &FF,
      const Eigen::MatrixBase<DerivedFFi> &FFi);

  template<typename DerivedF,
           typename DerivedFF,
           typename DerivedFFi>
  IGL_INLINE int triangle_tuple_get_face(
      const int &fi, const int &ei, const bool &along,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedFF> &FF,
      const Eigen::MatrixBase<DerivedFFi> &FFi);

  // Advanced Navigation
  /*!
     * Returns the next edge skipping the border
     *      _________
     *     /\ c | b /\
     *    /  \  |  /  \
     *   / d  \ | / a  \
     *  /______\|/______\
     *          v
     * In this example, if a and d are of-border and the pos is iterating counterclockwise, this method iterate through the faces incident on vertex v,
     * producing the sequence a, b, c, d, a, b, c, ...
     */
  template<typename DerivedF,
           typename DerivedFF,
           typename DerivedFFi>
  IGL_INLINE bool triangle_tuple_next_in_one_ring(
      int &fi, int &ei, bool &along,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedFF> &FF,
      const Eigen::MatrixBase<DerivedFFi> &FFi);

  // helper boolean functions
  template<typename DerivedF,
           typename DerivedFF,
           typename DerivedFFi>
  IGL_INLINE bool triangle_tuple_is_on_boundary(
      const int &fi, const int &ei, const bool &along,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedFF> &FF,
      const Eigen::MatrixBase<DerivedFFi> &FFi);

  IGL_INLINE bool triangle_tuples_equal(
      const int &f1, const int &e1, const bool &a1,
      const int &f2, const int &e2, const bool &a2);
}

#ifndef IGL_STATIC_LIBRARY
#  include "triangle_tuple.cpp"
#endif
#endif
