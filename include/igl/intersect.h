#ifndef IGL_INTERSECT_H
#define IGL_INTERSECT_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Determine the intersect between two sets of coefficients using ==
  // Templates:
  //   M  matrix type that implements indexing by global index M(i)
  // Inputs:
  //   A  matrix of coefficients
  //   B  matrix of coefficients
  // Output:
  //   C  matrix of elements appearing in both A and B, C is always resized to
  //   have a single column
  template <class M>
  IGL_INLINE void intersect(const M & A, const M & B, M & C);
  // Last argument as return
  template <class M>
  IGL_INLINE M intersect(const M & A, const M & B);
}
#ifdef IGL_HEADER_ONLY
#include "intersect.cpp"
#endif
#endif
