#ifndef IGL_AXIS_ANGLE_TO_QUAT_H
#define IGL_AXIS_ANGLE_TO_QUAT_H
#include "igl_inline.h"

#include "EPS.h"
#include <cmath>
namespace igl
{
  // Convert axis angle representation of a rotation to a quaternion
  // A Quaternion, q, is defined here as an arrays of four scalars (x,y,z,w),
  // such that q = x*i + y*j + z*k + w
  // Inputs:
  //   axis  3d vector
  //   angle  scalar
  // Outputs:
  //   quaternion
  template <typename Q_type>
  IGL_INLINE void axis_angle_to_quat(
    const Q_type *axis, 
    const Q_type angle,
    Q_type *out);
}

#ifdef IGL_HEADER_ONLY
#  include "axis_angle_to_quat.cpp"
#endif

#endif
