#ifndef MESH_BOOLEAN_CORK_H
#define MESH_BOOLEAN_CORK_H
#ifndef IGL_NO_CORK
#include "MeshBooleanType.h"
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <cork.h> // for consistent uint

namespace igl
{
  // Compute a boolean operation on two input meshes using the cork library.
  //
  // Inputs:
  //   VA  #VA by 3 list of vertex positions of first mesh
  //   FA  #FA by 3 list of triangle indices into VA
  //   VB  #VB by 3 list of vertex positions of second mesh
  //   FB  #FB by 3 list of triangle indices into VB
  //   type  of boolean operation see MeshBooleanType.h
  // Outputs:
  //   VC  #VC by 3 list of vertex positions of output mesh
  //   FC  #FC by 3 list of triangle indices into VC
  template <
    typename DerivedVA,
    typename DerivedFA,
    typename DerivedVB,
    typename DerivedFB,
    typename DerivedVC,
    typename DerivedFC>
  IGL_INLINE void mesh_boolean_cork(
    const Eigen::PlainObjectBase<DerivedVA > & VA,
    const Eigen::PlainObjectBase<DerivedFA > & FA,
    const Eigen::PlainObjectBase<DerivedVB > & VB,
    const Eigen::PlainObjectBase<DerivedFB > & FB,
    const MeshBooleanType & type,
    Eigen::PlainObjectBase<DerivedVC > & VC,
    Eigen::PlainObjectBase<DerivedFC > & FC);
}

#ifndef IGL_STATIC_LIBRARY
#  include "mesh_boolean_cork.cpp"
#endif
#endif

#endif
