#ifndef IGL_FROM_CORK_MESH_H
#define IGL_FROM_CORK_MESH_H
#ifndef IGL_NO_CORK
#include <igl/igl_inline.h>
#include <cork.h>
#include <Eigen/Core>
namespace igl
{
  // Convert cork's triangle mesh representation to a (V,F) mesh.
  //
  // Inputs:
  //   mesh  cork representation of mesh
  // Outputs:
  //   V  #V by 3 list of vertex positions
  //   F  #F by 3 list of triangle indices into V
  template <
    typename DerivedV,
    typename DerivedF>
  IGL_INLINE void from_cork_mesh(
    const CorkTriMesh & mesh,
    Eigen::PlainObjectBase<DerivedV > & V,
    Eigen::PlainObjectBase<DerivedF > & F);
}
#ifndef IGL_STATIC_LIBRARY
#  include "from_cork_mesh.cpp"
#endif
#endif
#endif
