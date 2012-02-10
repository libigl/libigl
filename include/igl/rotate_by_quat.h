#ifndef IGL_ROTATE_BY_QUAT_H
#define IGL_ROTATE_BY_QUAT_H
#include "igl_inline.h"

namespace igl
{
  // Compute rotation of a given vector/point by a quaternion
  // A Quaternion, q, is defined here as an arrays of four scalars (x,y,z,w),
  // such that q = x*i + y*j + z*k + w
  // Inputs:
  //   v  input 3d point/vector
  //   q  input quaternion
  // Outputs:
  //   out  result of rotation, allowed to be same as v
  template <typename Q_type>
  IGL_INLINE void rotate_by_quat(
    const Q_type *v,
    const Q_type *q, 
    Q_type *out);
};

#ifdef IGL_HEADER_ONLY
#  include "rotate_by_quat.cpp"
#endif

#endif
