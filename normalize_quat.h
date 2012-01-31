#ifndef IGL_NORMALIZE_QUAT_H
#define IGL_NORMALIZE_QUAT_H
#include "igl_inline.h"

namespace igl
{
  // Normalize a quaternion
  // A Quaternion, q, is defined here as an arrays of four scalars (x,y,z,w),
  // such that q = x*i + y*j + z*k + w
  // Inputs:
  //   q  input quaternion
  // Outputs:
  //   out  result of normalization, allowed to be same as q
  // Returns true on success, false if len(q) < EPS
  template <typename Q_type>
  IGL_INLINE bool normalize_quat(
    const Q_type *q,
    Q_type *out);
};

#ifdef IGL_HEADER_ONLY
#  include "normalize_quat.cpp"
#endif

#endif
