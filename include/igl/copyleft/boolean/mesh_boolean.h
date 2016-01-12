// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
//                    Qingnan Zhou <qnzhou@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
//
#ifndef IGL_COPYLEFT_BOOLEAN_MESH_BOOLEAN_H
#define IGL_COPYLEFT_BOOLEAN_MESH_BOOLEAN_H

#include "../../igl_inline.h"
#include "MeshBooleanType.h"
#include <Eigen/Core>
#include <functional>

namespace igl
{
  namespace copyleft
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
