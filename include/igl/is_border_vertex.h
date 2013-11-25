#ifndef IGL_IS_BORDER_VERTEX_H
#define IGL_IS_BORDER_VERTEX_H
#include "igl_inline.h"

#include <Eigen/Core>
#include <vector>

namespace igl 
{
  // Determine vertices on open boundary of a (manifold) mesh with triangle
  // faces F
  //
  // Inputs:
  //   V  #V by dim list of vertex positions 
  //   F  #F by 3 list of triangle indices
  // Returns vector of indices of vertices on open boundary of F
  //
  // Known Bugs: does not depend on V
  // 
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE std::vector<bool> is_border_vertex(
   const Eigen::PlainObjectBase<DerivedV> &V,
   const Eigen::PlainObjectBase<DerivedF> &F);
}

#ifdef IGL_HEADER_ONLY
#  include "is_border_vertex.cpp"
#endif

#endif
