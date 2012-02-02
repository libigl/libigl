#ifndef IGL_PER_FACE_NORMALS_H
#define IGL_PER_FACE_NORMALS_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Compute face normals via vertex position list, face list
  // Inputs:
  //   V  #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F  #F by 3 eigne Matrix of face (triangle) indices
  // Output:
  //   N  #F by 3 eigen Matrix of mesh face (triangle) 3D normals
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE void per_face_normals(
                                     const Eigen::PlainObjectBase<DerivedV>& V,
                                     const Eigen::PlainObjectBase<DerivedF>& F,
                                     Eigen::PlainObjectBase<DerivedV> & N);
}

#ifdef IGL_HEADER_ONLY
#  include "per_face_normals.cpp"
#endif

#endif
