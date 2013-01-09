#include "EPS.h"

template <> IGL_INLINE float igl::EPS()
{
  return igl::FLOAT_EPS;
}
template <> IGL_INLINE double igl::EPS()
{
  return igl::DOUBLE_EPS;
}

template <> IGL_INLINE float igl::EPS_SQ()
{
  return igl::FLOAT_EPS_SQ;
}
template <> IGL_INLINE double igl::EPS_SQ()
{
  return igl::DOUBLE_EPS_SQ;
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
