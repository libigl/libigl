// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
//                    Qingnan Zhou <qnzhou@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
//
#ifndef IGL_COPYLEFT_CGAL_MESH_BOOLEAN_H
#define IGL_COPYLEFT_CGAL_MESH_BOOLEAN_H

#include "../../igl_inline.h"
#include "../../MeshBooleanType.h"
#include <Eigen/Core>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <functional>

namespace igl
{
  namespace copyleft
  {
    namespace cgal
    {
      //  MESH_BOOLEAN Compute boolean csg operations on "solid", consistently
      //  oriented meshes.
      //
      //  Inputs:
      //    VA  #VA by 3 list of vertex positions of first mesh
      //    FA  #FA by 3 list of triangle indices into VA
      //    VB  #VB by 3 list of vertex positions of second mesh
      //    FB  #FB by 3 list of triangle indices into VB
      //    wind_num_op  function handle for filtering winding numbers from
      //      tuples of integer values to [0,1] outside/inside values
      //    keep  function handle for determining if a patch should be "kept"
      //      in the output based on the winding number on either side
      //    resolve_fun  function handle for computing resolve of a
      //      self-intersections of a mesh and outputting the new mesh.
      //  Outputs:
      //    VC  #VC by 3 list of vertex positions of boolean result mesh
      //    FC  #FC by 3 list of triangle indices into VC
      //    J  #FC list of indices into [FA;FB] revealing "birth" facet
      //  Returns true iff inputs induce a piecewise constant winding number
      //    field
      //
      //  See also: mesh_boolean_cork, intersect_other,
      //  remesh_self_intersections
      template <
        typename DerivedVA,
        typename DerivedFA,
        typename DerivedVB,
        typename DerivedFB,
        typename WindingNumberOp,
        typename KeepFunc,
        typename ResolveFunc,
        typename DerivedVC,
        typename DerivedFC,
        typename DerivedJ>
      IGL_INLINE bool mesh_boolean(
          const Eigen::PlainObjectBase<DerivedVA> & VA,
          const Eigen::PlainObjectBase<DerivedFA> & FA,
          const Eigen::PlainObjectBase<DerivedVB> & VB,
          const Eigen::PlainObjectBase<DerivedFB> & FB,
          const WindingNumberOp& wind_num_op,
          const KeepFunc& keep,
          const ResolveFunc& resolve_fun,
          Eigen::PlainObjectBase<DerivedVC > & VC,
          Eigen::PlainObjectBase<DerivedFC > & FC,
          Eigen::PlainObjectBase<DerivedJ > & J);
      ////  MESH_BOOLEAN Compute boolean csg operations on "solid", consistently
      ////  oriented meshes.
      ////
      ////  Inputs:
      ////    Vlist  k-long list of lists of mesh vertex positions
      ////    Flist  k-long list of lists of mesh face indices, so that Flist[i] indexes
      ////      vertices in Vlist[i]
      ////    wind_num_op  function handle for filtering winding numbers from
      ////      n-tuples of integer values to [0,1] outside/inside values
      ////    keep  function handle for determining if a patch should be "kept"
      ////      in the output based on the winding number on either side
      ////    resolve_fun  function handle for computing resolve of a
      ////      self-intersections of a mesh and outputting the new mesh.
      ////  Outputs:
      ////    VC  #VC by 3 list of vertex positions of boolean result mesh
      ////    FC  #FC by 3 list of triangle indices into VC
      ////    J  #FC list of indices into [Flist[0];Flist[1];...;Flist[k]]
      ////      revealing "birth" facet
      ////  Returns true iff inputs induce a piecewise constant winding number
      ////    field
      ////
      ////  See also: mesh_boolean_cork, intersect_other,
      ////  remesh_self_intersections
      //template <
      //  typename DerivedV,
      //  typename DerivedF,
      //  typename WindingNumberOp,
      //  typename KeepFunc,
      //  typename ResolveFunc,
      //  typename DerivedVC,
      //  typename DerivedFC,
      //  typename DerivedJ>
      //IGL_INLINE bool mesh_boolean(
      //    const std::vector<Eigen::PlainObjectBase<DerivedVA> > & Vlist,
      //    const std::vector<Eigen::PlainObjectBase<DerivedFA> > & Flist,
      //    const WindingNumberOp& wind_num_op,
      //    const KeepFunc& keep,
      //    const ResolveFunc& resolve_fun,
      //    Eigen::PlainObjectBase<DerivedVC > & VC,
      //    Eigen::PlainObjectBase<DerivedFC > & FC,
      //    Eigen::PlainObjectBase<DerivedJ > & J);
      // Given a resolved mesh (V,F), containing no self-intersections and a
      // list of birth parents.
      //
      // Inputs:
      //   V  #V by 3 list of resolved mesh vertex positions
      //   F  #F by 3 list of resolve mesh face indices 
      //   CJ #F list of birth parents 
      //   sizes  #inputs list of sizes so that sizes(i) is the #faces in the
      //     ith input
      //    wind_num_op  function handle for filtering winding numbers from
      //      tuples of integer values to [0,1] outside/inside values
      //    keep  function handle for determining if a patch should be "kept"
      //      in the output based on the winding number on either side
      //  Outputs:
      //    VC  #VC by 3 list of vertex positions of boolean result mesh
      //    FC  #FC by 3 list of triangle indices into VC
      //    J  #FC list of birth parent indices
      // 
      template <
        typename Derivedsizes,
        typename WindingNumberOp,
        typename KeepFunc,
        typename DerivedVC,
        typename DerivedFC,
        typename DerivedJ>
      IGL_INLINE bool mesh_boolean(
          const Eigen::Matrix<
            CGAL::Epeck::FT,
            Eigen::Dynamic,
            Eigen::Dynamic,
            DerivedVC::IsRowMajor> V,
          const Eigen::PlainObjectBase<DerivedFC > & F,
          const Eigen::Matrix<
            typename DerivedJ::Scalar,
            Eigen::Dynamic,
            1> & CJ,
          const Eigen::PlainObjectBase<Derivedsizes> & sizes,
          const WindingNumberOp& wind_num_op,
          const KeepFunc& keep,
          Eigen::PlainObjectBase<DerivedVC > & VC,
          Eigen::PlainObjectBase<DerivedFC > & FC,
          Eigen::PlainObjectBase<DerivedJ > & J);
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
      //  Returns true if inputs induce a piecewise constant winding number
      //    field and type is valid.
      //
      //  See also: mesh_boolean_cork, intersect_other,
      //  remesh_self_intersections
      template <
        typename DerivedVA,
        typename DerivedFA,
        typename DerivedVB,
        typename DerivedFB,
        typename ResolveFunc,
        typename DerivedVC,
        typename DerivedFC,
        typename DerivedJ>
      IGL_INLINE bool mesh_boolean(
          const Eigen::PlainObjectBase<DerivedVA > & VA,
          const Eigen::PlainObjectBase<DerivedFA > & FA,
          const Eigen::PlainObjectBase<DerivedVB > & VB,
          const Eigen::PlainObjectBase<DerivedFB > & FB,
          const MeshBooleanType & type,
        const ResolveFunc& resolve_func,
          Eigen::PlainObjectBase<DerivedVC > & VC,
          Eigen::PlainObjectBase<DerivedFC > & FC,
          Eigen::PlainObjectBase<DerivedJ > & J);

      //  Inputs:
      //    VA  #VA by 3 list of vertex positions of first mesh
      //    FA  #FA by 3 list of triangle indices into VA
      //    VB  #VB by 3 list of vertex positions of second mesh
      //    FB  #FB by 3 list of triangle indices into VB
      //    type  type of boolean operation
      //  Outputs:
      //    VC  #VC by 3 list of vertex positions of boolean result mesh
      //    FC  #FC by 3 list of triangle indices into VC
      //    J  #FC list of indices into [FA;FB] revealing "birth" facet
      //  Returns true if inputs induce a piecewise constant winding number
      //  field and type is valid
      //
      //  See also: mesh_boolean_cork, intersect_other,
      //  remesh_self_intersections
      template <
        typename DerivedVA,
        typename DerivedFA,
        typename DerivedVB,
        typename DerivedFB,
        typename DerivedVC,
        typename DerivedFC,
        typename DerivedJ>
      IGL_INLINE bool mesh_boolean(
        const Eigen::PlainObjectBase<DerivedVA > & VA,
        const Eigen::PlainObjectBase<DerivedFA > & FA,
        const Eigen::PlainObjectBase<DerivedVB > & VB,
        const Eigen::PlainObjectBase<DerivedFB > & FB,
        const MeshBooleanType & type,
        Eigen::PlainObjectBase<DerivedVC > & VC,
        Eigen::PlainObjectBase<DerivedFC > & FC,
        Eigen::PlainObjectBase<DerivedJ > & J);

      //  Inputs:
      //    VA  #VA by 3 list of vertex positions of first mesh
      //    FA  #FA by 3 list of triangle indices into VA
      //    VB  #VB by 3 list of vertex positions of second mesh
      //    FB  #FB by 3 list of triangle indices into VB
      //    type  type of boolean operation
      //  Outputs:
      //    VC  #VC by 3 list of vertex positions of boolean result mesh
      //    FC  #FC by 3 list of triangle indices into VC
      //  Returns true ff inputs induce a piecewise constant winding number
      //    field and type is valid
      template <
        typename DerivedVA,
        typename DerivedFA,
        typename DerivedVB,
        typename DerivedFB,
        typename DerivedVC,
        typename DerivedFC>
      IGL_INLINE bool mesh_boolean(
          const Eigen::PlainObjectBase<DerivedVA > & VA,
          const Eigen::PlainObjectBase<DerivedFA > & FA,
          const Eigen::PlainObjectBase<DerivedVB > & VB,
          const Eigen::PlainObjectBase<DerivedFB > & FB,
          const MeshBooleanType & type,
          Eigen::PlainObjectBase<DerivedVC > & VC,
          Eigen::PlainObjectBase<DerivedFC > & FC);
    }
  }
}

#ifndef IGL_STATIC_LIBRARY
#  include "mesh_boolean.cpp"
#endif

#endif
