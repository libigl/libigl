#ifndef IGL_TO_CORK_MESH_H
#define IGL_TO_CORK_MESH_H
#ifndef IGL_NO_CORK
#include <igl/igl_inline.h>
#include <cork.h>
#include <Eigen/Core>
namespace igl
{
  // Convert a (V,F) mesh to a cork's triangle mesh representation.
  //
  // Inputs:
  //   V  #V by 3 list of vertex positions
  //   F  #F by 3 list of triangle indices into V
  // Outputs:
  //   mesh  cork representation of mesh
  template <
    typename DerivedV,
    typename DerivedF>
  IGL_INLINE void to_cork_mesh(
    const Eigen::PlainObjectBase<DerivedV > & V,
    const Eigen::PlainObjectBase<DerivedF > & F,
    CorkTriMesh & mesh);
}
#ifndef IGL_STATIC_LIBRARY
#  include "to_cork_mesh.cpp"
#endif
#endif
#endif
