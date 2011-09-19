#ifndef IGL_ROTATE_BY_QUAT_H
#define IGL_ROTATE_BY_QUAT_H

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
  inline void rotate_by_quat(
    const Q_type *v,
    const Q_type *q, 
    Q_type *out);
};

// Implementation
#include "quat_conjugate.h"
#include "quat_mult.h"
#include "normalize_quat.h"
#include <cassert>

template <typename Q_type>
inline void igl::rotate_by_quat(
  const Q_type *v,
  const Q_type *q,
  Q_type *out)
{
  // Quaternion form of v, copy data in v, (as a result out can be same pointer
  // as v)
  Q_type quat_v[4] = {v[0],v[1],v[2],0};

  // normalize input 
  Q_type normalized_q[4];
  bool normalized = igl::normalize_quat<Q_type>(q,normalized_q);
  assert(normalized);

  // Conjugate of q
  Q_type q_conj[4];
  igl::quat_conjugate<Q_type>(normalized_q,q_conj);

  // Rotate of vector v by quaternion q is:
  // q*v*conj(q)
  // Compute q*v
  Q_type q_mult_quat_v[4];
  igl::quat_mult<Q_type>(normalized_q,quat_v,q_mult_quat_v);
  // Compute (q*v) * conj(q)
  igl::quat_mult<Q_type>(q_mult_quat_v,q_conj,out);
}

#endif

