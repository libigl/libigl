#ifndef IGL_EPS_H
#define IGL_EPS_H
#include "igl_inline.h"
// Define a standard value for double epsilon
namespace igl
{
  const double DOUBLE_EPS    = 1.0e-14;
  const double DOUBLE_EPS_SQ = 1.0e-28;
  const float FLOAT_EPS    = 1.0e-7;
  const float FLOAT_EPS_SQ = 1.0e-14;
  // Function returning EPS for corresponding type
  template <typename S_type> IGL_INLINE S_type EPS();
  template <typename S_type> IGL_INLINE S_type EPS_SQ();
  // Template specializations for float and double
  template <> IGL_INLINE float EPS<float>();
  template <> IGL_INLINE double EPS<double>();
  template <> IGL_INLINE float EPS_SQ<float>();
  template <> IGL_INLINE double EPS_SQ<double>();
}

#ifdef IGL_HEADER_ONLY
#  include "EPS.cpp"
#endif

#endif
