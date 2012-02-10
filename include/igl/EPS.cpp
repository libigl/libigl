#include "EPS.h"

template <> IGL_INLINE float igl::EPS()
{
  return igl::FLOAT_EPS;
}

template <> IGL_INLINE double igl::EPS()
{
  return igl::DOUBLE_EPS;
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
