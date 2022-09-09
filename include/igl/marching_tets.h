// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Francis Williams <francis@fwilliams.info>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_MARCHING_TETS_H
#define IGL_MARCHING_TETS_H

#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace igl {
  // marching_tets( TV, TT, S, isovalue, SV, SF, J, BC)
  //
  // performs the marching tetrahedra algorithm on a tet mesh defined by TV and
  // TT with scalar values defined at each vertex in TV. The output is a
  // triangle mesh approximating the isosurface coresponding to the value
  // isovalue.
  //
  // Input:
  //  TV  #tet_vertices x 3 array -- The vertices of the tetrahedral mesh
  //  TT  #tets x 4 array --  The indexes of each tet in the tetrahedral mesh
  //  S  #tet_vertices x 1 array -- The values defined on each tet vertex
  //  isovalue  scalar -- The isovalue of the level set we want to compute
  //
  // Output:
  //  SV  #SV x 3 array -- The vertices of the output level surface mesh
  //  SF  #SF x 3 array -- The face indexes of the output level surface mesh
  //  J   #SF list of indices into TT revealing which tet each face comes from
  //  BC  #SV x #TV list of barycentric coordinates so that SV = BC*TV
  template <typename DerivedTV,
            typename DerivedTT,
            typename DerivedS,
            typename DerivedSV,
            typename DerivedSF,
            typename DerivedJ,
            typename BCType>
  IGL_INLINE void marching_tets(
      const Eigen::MatrixBase<DerivedTV>& TV,
      const Eigen::MatrixBase<DerivedTT>& TT,
      const Eigen::MatrixBase<DerivedS>& S,
      double isovalue,
      Eigen::PlainObjectBase<DerivedSV>& SV,
      Eigen::PlainObjectBase<DerivedSF>& SF,
      Eigen::PlainObjectBase<DerivedJ>& J,
      Eigen::SparseMatrix<BCType>& BC);

  // marching_tets( TV, TT, S, SV, SF, J, BC)
  //
  // Performs the marching tetrahedra algorithm on a tet mesh defined by TV and
  // TT with scalar values defined at each vertex in TV. The output is a
  // triangle mesh approximating the isosurface coresponding to an isovalue of 0.
  //
  // Input:
  //  TV  #tet_vertices x 3 array -- The vertices of the tetrahedral mesh
  //  TT  #tets x 4 array --  The indexes of each tet in the tetrahedral mesh
  //  S  #tet_vertices x 1 array -- The values defined on each tet vertex
  //  isovalue  scalar -- The isovalue of the level set we want to compute
  //
  // Output:
  //  SV  #SV x 3 array -- The vertices of the output level surface mesh
  //  SF  #SF x 3 array -- The face indexes of the output level surface mesh
  //  J   #SF list of indices into TT revealing which tet each face comes from
  //  BC  #SV x #TV list of barycentric coordinates so that SV = BC*TV
  template <typename DerivedTV,
            typename DerivedTT,
            typename DerivedS,
            typename DerivedSV,
            typename DerivedSF,
            typename DerivedJ,
            typename BCType>
  IGL_INLINE void marching_tets(
      const Eigen::MatrixBase<DerivedTV>& TV,
      const Eigen::MatrixBase<DerivedTT>& TT,
      const Eigen::MatrixBase<DerivedS>& S,
      Eigen::PlainObjectBase<DerivedSV>& SV,
      Eigen::PlainObjectBase<DerivedSF>& SF,
      Eigen::PlainObjectBase<DerivedJ>& J,
      Eigen::SparseMatrix<BCType>& BC) {
    return igl::marching_tets(TV, TT, S, 0.0, SV, SF, J, BC);
  }

  // marching_tets( TV, TT, S, isovalue, SV, SF, J)
  //
  // performs the marching tetrahedra algorithm on a tet mesh defined by TV and
  // TT with scalar values defined at each vertex in TV. The output is a
  // triangle mesh approximating the isosurface coresponding to the value
  // isovalue.
  //
  // Input:
  //  TV  #tet_vertices x 3 array -- The vertices of the tetrahedral mesh
  //  TT  #tets x 4 array --  The indexes of each tet in the tetrahedral mesh
  //  S  #tet_vertices x 1 array -- The values defined on each tet vertex
  //  isovalue  scalar -- The isovalue of the level set we want to compute
  //
  // Output:
  //  SV  #SV x 3 array -- The vertices of the output level surface mesh
  //  SF  #SF x 3 array -- The face indexes of the output level surface mesh
  //  J   #SF list of indices into TT revealing which tet each face comes from
  template <typename DerivedTV,
            typename DerivedTT,
            typename DerivedS,
            typename DerivedSV,
            typename DerivedSF,
            typename DerivedJ>
  IGL_INLINE void marching_tets(
      const Eigen::MatrixBase<DerivedTV>& TV,
      const Eigen::MatrixBase<DerivedTT>& TT,
      const Eigen::MatrixBase<DerivedS>& S,
      double isovalue,
      Eigen::PlainObjectBase<DerivedSV>& SV,
      Eigen::PlainObjectBase<DerivedSF>& SF,
      Eigen::PlainObjectBase<DerivedJ>& J) {
    Eigen::SparseMatrix<double> _BC;
    return igl::marching_tets(TV, TT, S, isovalue, SV, SF, J, _BC);
  }

  // marching_tets( TV, TT, S, isovalue, SV, SF, BC)
  //
  // performs the marching tetrahedra algorithm on a tet mesh defined by TV and
  // TT with scalar values defined at each vertex in TV. The output is a
  // triangle mesh approximating the isosurface coresponding to the value
  // isovalue.
  //
  // Input:
  //  TV  #tet_vertices x 3 array -- The vertices of the tetrahedral mesh
  //  TT  #tets x 4 array --  The indexes of each tet in the tetrahedral mesh
  //  S  #tet_vertices x 1 array -- The values defined on each tet vertex
  //  isovalue  scalar -- The isovalue of the level set we want to compute
  //
  // Output:
  //  SV  #SV x 3 array -- The vertices of the output level surface mesh
  //  SF  #SF x 3 array -- The face indexes of the output level surface mesh
  //  BC  #SV x #TV list of barycentric coordinates so that SV = BC*TV
  template <typename DerivedTV,
            typename DerivedTT,
            typename DerivedS,
            typename DerivedSV,
            typename DerivedSF,
            typename BCType>
  IGL_INLINE void marching_tets(
      const Eigen::MatrixBase<DerivedTV>& TV,
      const Eigen::MatrixBase<DerivedTT>& TT,
      const Eigen::MatrixBase<DerivedS>& S,
      double isovalue,
      Eigen::PlainObjectBase<DerivedSV>& SV,
      Eigen::PlainObjectBase<DerivedSF>& SF,
      Eigen::SparseMatrix<BCType>& BC) {
    Eigen::VectorXi _J;
    return igl::marching_tets(TV, TT, S, isovalue, SV, SF, _J, BC);
  }

  // marching_tets( TV, TT, S, isovalue, SV, SF)
  //
  // performs the marching tetrahedra algorithm on a tet mesh defined by TV and
  // TT with scalar values defined at each vertex in TV. The output is a
  // triangle mesh approximating the isosurface coresponding to the value
  // isovalue.
  //
  // Input:
  //  TV  #tet_vertices x 3 array -- The vertices of the tetrahedral mesh
  //  TT  #tets x 4 array --  The indexes of each tet in the tetrahedral mesh
  //  S  #tet_vertices x 1 array -- The values defined on each tet vertex
  //  isovalue  scalar -- The isovalue of the level set we want to compute
  //
  // Output:
  //  SV  #SV x 3 array -- The vertices of the output level surface mesh
  //  SF  #SF x 3 array -- The face indexes of the output level surface mesh
  template <typename DerivedTV,
            typename DerivedTT,
            typename DerivedS,
            typename DerivedSV,
            typename DerivedSF>
  IGL_INLINE void marching_tets(
      const Eigen::MatrixBase<DerivedTV>& TV,
      const Eigen::MatrixBase<DerivedTT>& TT,
      const Eigen::MatrixBase<DerivedS>& S,
      double isovalue,
      Eigen::PlainObjectBase<DerivedSV>& SV,
      Eigen::PlainObjectBase<DerivedSF>& SF) {
    Eigen::VectorXi _J;
    Eigen::SparseMatrix<double> _BC;
    return igl::marching_tets(TV, TT, S, isovalue, SV, SF, _J, _BC);
  }

}

#ifndef IGL_STATIC_LIBRARY
#  include "marching_tets.cpp"
#endif

#endif // IGL_MARCHING_TETS_H
