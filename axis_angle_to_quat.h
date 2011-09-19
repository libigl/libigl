#ifndef IGL_AXIS_ANGLE_TO_QUAT_H
#define IGL_AXIS_ANGLE_TO_QUAT_H

#include <EPS.h>
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
  inline void axis_angle_to_quat(
    const Q_type *axis, 
    const Q_type angle,
    Q_type *out);
}

// Implementation
// http://www.antisphere.com/Wiki/tools:anttweakbar
template <typename Q_type>
inline void igl::axis_angle_to_quat(
  const Q_type *axis, 
  const Q_type angle,
  Q_type *out)
{
    Q_type n = axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2];
    if( fabs(n)>igl::EPS<Q_type>())
    {
        Q_type f = 0.5*angle;
        out[3] = cos(f);
        f = sin(f)/sqrt(n);
        out[0] = axis[0]*f;
        out[1] = axis[1]*f;
        out[2] = axis[2]*f;
    }
    else
    {
        out[3] = 1.0;
        out[0] = out[1] = out[2] = 0.0;
    }
}

#endif
