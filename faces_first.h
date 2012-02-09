#ifndef IGL_FACES_FIRST_H
#define IGL_FACES_FIRST_H
#include "igl_inline.h"
namespace igl
{
  // FACES_FIRST Reorder vertices so that vertices in face list come before
  // vertices that don't appear in the face list. This is especially useful if
  // the face list contains only surface faces and you want surface vertices
  // listed before internal vertices
  //
  // [RV,RT,RF,IM] = faces_first(V,T,F);
  //
  // Templates:
  //   MatV  matrix for vertex positions, e.g. MatrixXd
  //   MatF  matrix for face indices, e.g. MatrixXi
  //   VecI  vector for index map, e.g. VectorXi
  // Input:
  //  V  # vertices by 3 vertex positions
  //  F  # faces by 3 list of face indices
  // Output: 
  //  RV  # vertices by 3 vertex positions, order such that if the jth vertex is
  //    some face in F, and the kth vertex is not then j comes before k
  //  RF  # faces by 3 list of face indices, reindexed to use RV
  //  IM  # faces by 1 list of indices such that: RF = IM(F) and RT = IM(T)
  //    and RV(IM,:) = V
  //
  template <typename MatV, typename MatF, typename VecI>
  IGL_INLINE void faces_first(
    const MatV & V, 
    const MatF & F, 
    MatV & RV, 
    MatF & RF, 
    VecI & IM);
  // Virtual "in place" wrapper
  template <typename MatV, typename MatF, typename VecI>
  IGL_INLINE void faces_first(
    MatV & V, 
    MatF & F, 
    VecI & IM);
}

#ifdef IGL_HEADER_ONLY
#  include "faces_first.cpp"
#endif

#endif
