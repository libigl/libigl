#ifndef IGL_QUAT_TO_AXIS_ANGLE_H
#define IGL_QUAT_TO_AXIS_ANGLE_H
#include "igl_inline.h"

namespace igl
{
  // Convert quat representation of a rotation to axis angle
  // A Quaternion, q, is defined here as an arrays of four scalars (x,y,z,w),
  // such that q = x*i + y*j + z*k + w
  // Inputs:
  //   q quaternion
  // Outputs:
  //   axis  3d vector
  //   angle  scalar in radians
  template <typename Q_type>
  IGL_INLINE void quat_to_axis_angle(
    const Q_type *q,
    Q_type *axis, 
    Q_type & angle);
  // Wrapper with angle in degrees
  template <typename Q_type>
  IGL_INLINE void quat_to_axis_angle_deg(
    const Q_type *q,
    Q_type *axis, 
    Q_type & angle);
}

#ifdef IGL_HEADER_ONLY
#  include "quat_to_axis_angle.cpp"
#endif

#endif

