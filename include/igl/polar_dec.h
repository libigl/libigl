#ifndef IGL_POLAR_DEC
#define IGL_POLAR_DEC
#include "igl_inline.h"

namespace igl
{
  // Computes the polar decomposition (R,T) of a matrix A
  // Inputs:
  //   A  3 by 3 matrix to be decomposed
  // Outputs:
  //   R  3 by 3 rotation matrix part of decomposition
  //   T  3 by 3 stretch matrix part of decomposition
  //
  // Note: I'm not sure if this implementation is check against reflections in R
  // Note: It is not
  //
  template<typename Mat>
  IGL_INLINE void polar_dec(const Mat& A, Mat& R, Mat& T);
}
#ifdef IGL_HEADER_ONLY
#  include "polar_dec.cpp"
#endif
#endif

