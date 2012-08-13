#ifndef IGL_PER_VERTEX_ATTRIBUTE_SMOOTHING_H
#define IGL_PER_VERTEX_ATTRIBUTE_SMOOTHING_H
#include "igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  // Smooth vertex attributes using uniform Laplacian
  // Inputs:
  //   Ain  #V by #A eigen Matrix of mesh vertex attributes (each vertex has #A attributes)
  //   F    #F by 3 eigne Matrix of face (triangle) indices
  // Output:
  //   Aout #V by #A eigen Matrix of mesh vertex attributes
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE void per_vertex_attribute_smoothing(
    const Eigen::PlainObjectBase<DerivedV>& Ain,
    const Eigen::PlainObjectBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedV> & Aout);
}

#ifdef IGL_HEADER_ONLY
#  include "per_vertex_attribute_smoothing.cpp"
#endif

#endif
