#include "canonical_quaternions.h"

template <> IGL_INLINE const float igl::CANONICAL_VIEW_QUAT<float>(int i, int j)
{
  return (float)igl::CANONICAL_VIEW_QUAT_F[i][j];
}
template <> IGL_INLINE const double igl::CANONICAL_VIEW_QUAT<double>(int i, int j)
{
  return (double)igl::CANONICAL_VIEW_QUAT_D[i][j];
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
