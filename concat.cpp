#include "concat.h"

#include <cstdio>

template <typename T>
IGL_INLINE void igl::concat(
                   const T A, 
                   const T B,
                   const bool horiz,                 
                   T& O)
{
  if (horiz)
  {
    // O = [A,B]
    assert(A.rows() == B.rows());
    O = T(A.rows(),A.cols()+B.cols());
    O << A,B;
  }
  else
  {
    // O = [A;B]
    assert(A.cols() == B.cols());
    O = T(A.rows()+B.rows(),A.cols());
    O << A,B;
  }
}

template <typename T>
IGL_INLINE T igl::concat(
                const T A, 
                const T B,
                bool horiz
                )
{
  T O = T(1,1);
  concat(A,B,horiz,O);
  return O;
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
