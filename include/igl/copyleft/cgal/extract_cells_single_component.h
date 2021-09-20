// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Qingnan Zhou <qnzhou@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
//
#ifndef IGL_COPYLEFT_CGAL_EXTRACT_CELLS_SINGLE_COMPONENT_H
#define IGL_COPYLEFT_CGAL_EXTRACT_CELLS_SINGLE_COMPONENT_H

#include "../../igl_inline.h"
#include <Eigen/Core>
#include <vector>

namespace igl {
  namespace copyleft
  {
    namespace cgal
    {
      // Extract connected 3D space partitioned by mesh (V,F) composed of
      // **possibly multiple components** (the name of this function is
      // dubious).
      //
      // Inputs:
      //   V  #V by 3 array of vertices.
      //   F  #F by 3 array of faces.
      //   P  #F list of patch indices.
      //   E  #E by 2 array of vertex indices, one edge per row.
      //  uE    #uE by 2 list of vertex_indices, represents undirected edges.
      //  EMAP  #F*3 list of indices into uE.
      //  uEC  #uE+1 list of cumsums of directed edges sharing each unique edge
      //  uEE  #E list of indices into E (see `igl::unique_edge_map`)
      // Output:
      //  cells  #P by 2 array of cell indices.  cells(i,0) represents the
      //    cell index on the positive side of patch i, and cells(i,1)
      //    represents cell index of the negative side.
      template<
        typename DerivedV,
        typename DerivedF,
        typename DerivedP,
        typename DeriveduE,
        typename DerivedEMAP,
        typename DeriveduEC,
        typename DeriveduEE,
        typename DerivedC >
      IGL_INLINE size_t extract_cells_single_component(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        const Eigen::PlainObjectBase<DerivedP>& P,
        const Eigen::PlainObjectBase<DeriveduE>& uE,
        const Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
        const Eigen::PlainObjectBase<DeriveduEC>& uEC,
        const Eigen::PlainObjectBase<DeriveduEE>& uEE,
        Eigen::PlainObjectBase<DerivedC>& cells);
    }
  }
}

#ifndef IGL_STATIC_LIBRARY
#  include "extract_cells_single_component.cpp"
#endif
#endif

