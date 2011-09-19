#ifndef IGL_NORMALIZE_QUAT_H
#define IGL_NORMALIZE_QUAT_H

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
  inline bool normalize_quat(
    const Q_type *q,
    Q_type *out);
};

// Implementation
#include "EPS.h"

template <typename Q_type>
inline bool igl::normalize_quat(
  const Q_type *q,
  Q_type *out)
{
  // Get length
  Q_type len = sqrt(
    q[0]*q[0]+
    q[1]*q[1]+
    q[2]*q[2]+
    q[3]*q[3]);

  // Noramlize each coordinate
  out[0] = q[0]/len;
  out[1] = q[1]/len;
  out[2] = q[2]/len;
  out[3] = q[3]/len;

  // Test whether length was below Epsilon
  return (len > igl::EPS<Q_type>());
}

#endif


