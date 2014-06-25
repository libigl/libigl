#ifndef IGL_QUAD_PLANARITY_H
#define IGL_QUAD_PLANARITY_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Compute planarity of the faces of a quad mesh
  // Inputs:
  //   V  #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F  #F by 4 eigen Matrix of face (quad) indices
  // Output:
  //   P  #F by 1 eigen Matrix of mesh face (quad) planarities
  //
  template <typename DerivedV, typename DerivedF, typename DerivedP>
  IGL_INLINE void quad_planarity(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedP> & P);
}

#ifndef IGL_STATIC_LIBRARY
#  include "quad_planarity.cpp"
#endif

#endif
