#ifndef IGL_COTANGENT_H
#define IGL_COTANGENT_H
#include "igl_inline.h"
namespace igl
{
  // COTANGENT compute the cotangents of each angle in mesh (V,F)
  // 
  // Templates:
  //   MatV  vertex position matrix, e.g. Eigen::MatrixXd
  //   MatF  face index matrix, e.g. Eigen::MatrixXd
  //   MatC  cotangent weights matrix, e.g. Eigen::MatrixXd
  // Inputs:
  //   V  #V by dim list of rest domain positions
  //   F  #F by {3|4} list of {triangle|tetrahedra} indices into V
  // Outputs:
  //   C  #F by {3|6} list of cotangents corresponding angles
  //     for triangles, columns correspond to edges 23,31,12
  //     for tets, columns correspond to edges 23,31,12,41,42,43
  template <class MatV, class MatF, class MatC>
  IGL_INLINE void cotangent(const MatV & V, const MatF & F, MatC & C);
}

#ifdef IGL_HEADER_ONLY
#  include "cotangent.cpp"
#endif

#endif
