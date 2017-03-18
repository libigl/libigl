// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Zhongshi Jiang <zhongshi@cims.nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_TETRAHEDRON_TETRAHEDRON_ADJACENCY_H
#define IGL_TETRAHEDRON_TETRAHEDRON_ADJACENCY_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>

namespace igl
{
  // Constructs the tetrahedron-tetrahedron adjacency matrix for a given
  // tet mesh (V,T).
  //
  // Inputs:
  //   T  #T by 4 list of mesh faces (must be tetrahedrons)
  // Outputs:
  //   TT   #T by #4  adjacent matrix, the element i,j is the id of
  //          the tetrahedron adjacent to the j face of tetrahedron i
  //   TTif #T by #4  adjacent matrix, the element i,j is the id of
  //          face of the tetrahedron TT(i,j) adjacent with tetrahedron i
  //   TTie #T by #12 adjacent matrix, the element i, 3*j+k is the id of
  //          edge of the face TTif(i,j) of tet TT(i,j)
  //          adjacent with tetrahedron i
  // NOTE:
  //   Faces are indexed by opposite corner vertex and
  //   orientation is determined looking from the opposite vertex,
  //   For tet T = [v0,v1,v2,v3]:
  //    F[0] = v3,v2,v1;
  //    F[1] = v2,v3,v0;
  //    F[2] = v1,v0,v3;
  //    F[3] = v0,v1,v2;
  //   However, to be compatible with triangle_triangle_adjacency,
  //   edges are indexed with its starting vertex:
  //   the first edge for [v0,v1,v2] is [v0,v1] the second [v1,v2]
  //   and the third [v2,v3].
  //       this convention is DIFFERENT from cotmatrix_entries.h
  template<typename DerivedT, typename DerivedTT,
           typename DerivedTTif, typename DerivedTTie>
  IGL_INLINE void tetrahedron_tetrahedron_adjacency(
      const Eigen::MatrixBase<DerivedT> &T,
      Eigen::MatrixBase<DerivedTT> &TT,
      Eigen::MatrixBase<DerivedTTif> &TTif,
      Eigen::MatrixBase<DerivedTTie> &TTie);

  template<typename DerivedT, typename DerivedTT>
  IGL_INLINE void tetrahedron_tetrahedron_adjacency(
      const Eigen::MatrixBase<DerivedT> &T,
      Eigen::MatrixBase<DerivedTT> &TT);

  // Helper Matrix
  // Local Triangle-Triangle-Adjacency.
  // Face indexing matrix from the opposite corner
  // the matrix has the property that FFi is always identity with edge.
  // can be rephrased as (i-j+3+(i%2)*(2*j-2))%4
  static Eigen::Matrix<int, 4, 3> tetrahedron_local_FF =
      (Eigen::Matrix<int, 4, 3>() <<
          3, 2, 1,
          2, 3, 0,
          1, 0, 3,
          0, 1, 2).finished();
}

#ifndef IGL_STATIC_LIBRARY
#  include "tetrahedron_tetrahedron_adjacency.cpp"
#endif

#endif
