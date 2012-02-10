#ifndef IGL_SNAP_TO_CANONICAL_VIEW_QUAT_H
#define IGL_SNAP_TO_CANONICAL_VIEW_QUAT_H
#include "igl_inline.h"
// A Quaternion, q, is defined here as an arrays of four scalars (x,y,z,w),
// such that q = x*i + y*j + z*k + w
namespace igl
{
  // Snap a given quaternion to the "canonical quaternion" rotations.
  // Inputs:
  //   q  input quaternion
  //   threshold  threshold between 0 and 1, where 0 means
  template <typename Q_type>
  IGL_INLINE bool snap_to_canonical_view_quat(
    const Q_type q[4],
    const Q_type threshold,
    Q_type s[4]);
}

#ifdef IGL_HEADER_ONLY
#  include "snap_to_canonical_view_quat.cpp"
#endif

#endif
