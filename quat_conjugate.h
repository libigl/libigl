#ifndef IGL_QUAT_CONJUGATE_H
#define IGL_QUAT_CONJUGATE_H

namespace igl
{
  // Compute conjugate of given quaternion
  // http://en.wikipedia.org/wiki/Quaternion#Conjugation.2C_the_norm.2C_and_reciprocal
  // A Quaternion, q, is defined here as an arrays of four scalars (x,y,z,w),
  // such that q = x*i + y*j + z*k + w
  // Inputs:
  //   q1  input quaternion
  // Outputs:
  //   out  result of conjugation, allowed to be same as input
  template <typename Q_type>
  inline void quat_conjugate(
    const Q_type *q1, 
    Q_type *out);
};

// Implementation
template <typename Q_type>
inline void igl::quat_conjugate(
  const Q_type *q1, 
  Q_type *out)
{
  out[0] = -q1[0];
  out[1] = -q1[1];
  out[2] = -q1[2];
  out[3] = q1[3];
}

#endif

