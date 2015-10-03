// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_BOOLEAN_MESH_BOOLEAN_H
#define IGL_BOOLEAN_MESH_BOOLEAN_H

#include <igl/igl_inline.h>
#include "MeshBooleanType.h"
#include <Eigen/Core>
#include <functional>

namespace igl
{
  namespace boolean
  {
    //  MESH_BOOLEAN Compute boolean csg operations on "solid", consistently
    //  oriented meshes.
    // 
    //  Inputs:
    //    VA  #VA by 3 list of vertex positions of first mesh
    //    FA  #FA by 3 list of triangle indices into VA
    //    VB  #VB by 3 list of vertex positions of second mesh
    //    FB  #FB by 3 list of triangle indices into VB
    //    type  type of boolean operation
    //    resolve_fun  function handle for computing resolve of a
    //      self-intersections of a mesh and outputting the new mesh.
    //  Outputs:
    //    VC  #VC by 3 list of vertex positions of boolean result mesh
    //    FC  #FC by 3 list of triangle indices into VC
    //    J  #FC list of indices into [FA;FB] revealing "birth" facet
    //    I  #VA+#VB list of indices into SV (SV=[VA;VB;SVA;SVB], where SVA and
    //      SVB are the new vertices from resolving intersections) revealing
    //      "birth" vertices.
    //
    //  See also: mesh_boolean_cork, intersect_other, remesh_self_intersections
    //     
    template <
      typename DerivedVA,
      typename DerivedFA,
      typename DerivedVB,
      typename DerivedFB,
      typename DerivedVC,
      typename DerivedFC,
      typename DerivedJ,
      typename DerivedI>
    IGL_INLINE void mesh_boolean(
      const Eigen::PlainObjectBase<DerivedVA > & VA,
      const Eigen::PlainObjectBase<DerivedFA > & FA,
      const Eigen::PlainObjectBase<DerivedVB > & VB,
      const Eigen::PlainObjectBase<DerivedFB > & FB,
      const MeshBooleanType & type,
      const std::function<void(
        const Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3> &,
        const Eigen::Matrix<typename DerivedFC::Scalar,Eigen::Dynamic,3> &,
        Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3> &,
        Eigen::Matrix<typename DerivedFC::Scalar,Eigen::Dynamic,3> &,
        Eigen::Matrix<typename DerivedJ::Scalar,Eigen::Dynamic,1>&,
        Eigen::Matrix<typename DerivedI::Scalar,Eigen::Dynamic,1>&)> &
        resolve_fun,
      Eigen::PlainObjectBase<DerivedVC > & VC,
      Eigen::PlainObjectBase<DerivedFC > & FC,
      Eigen::PlainObjectBase<DerivedJ > & J,
      Eigen::PlainObjectBase<DerivedI > & I);
    template <
      typename DerivedVA,
      typename DerivedFA,
      typename DerivedVB,
      typename DerivedFB,
      typename DerivedVC,
      typename DerivedFC,
      typename DerivedJ,
      typename DerivedI>
    IGL_INLINE void mesh_boolean(
      const Eigen::PlainObjectBase<DerivedVA > & VA,
      const Eigen::PlainObjectBase<DerivedFA > & FA,
      const Eigen::PlainObjectBase<DerivedVB > & VB,
      const Eigen::PlainObjectBase<DerivedFB > & FB,
      const MeshBooleanType & type,
      Eigen::PlainObjectBase<DerivedVC > & VC,
      Eigen::PlainObjectBase<DerivedFC > & FC,
      Eigen::PlainObjectBase<DerivedJ > & J,
      Eigen::PlainObjectBase<DerivedI > & I);
    template <
      typename DerivedVA,
      typename DerivedFA,
      typename DerivedVB,
      typename DerivedFB,
      typename DerivedVC,
      typename DerivedFC>
    IGL_INLINE void mesh_boolean(
      const Eigen::PlainObjectBase<DerivedVA > & VA,
      const Eigen::PlainObjectBase<DerivedFA > & FA,
      const Eigen::PlainObjectBase<DerivedVB > & VB,
      const Eigen::PlainObjectBase<DerivedFB > & FB,
      const MeshBooleanType & type,
      Eigen::PlainObjectBase<DerivedVC > & VC,
      Eigen::PlainObjectBase<DerivedFC > & FC);

  }
}

#ifndef IGL_STATIC_LIBRARY
#  include "mesh_boolean.cpp"
#endif

#endif
