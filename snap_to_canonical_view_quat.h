#ifndef IGL_SNAP_TO_CANONICAL_QUAT_H
#define IGL_SNAP_TO_CANONICAL_QUAT_H
// A Quaternion, q, is defined here as an arrays of four scalars (x,y,z,w),
// such that q = x*i + y*j + z*k + w
namespace igl
{
  // Snap a given quaternion to the "canonical quaternion" rotations.
  // Inputs:
  //   q  input quaternion
  //   threshold  threshold between 0 and 1, where 0 means
  template <typename Q_type>
  inline bool snap_to_canonical_view_quat(
    const Q_type q[4],
    const Q_type threshold,
    Q_type s[4]);
}


// Implementation
#include "canonical_quaternions.h"
#include "normalize_quat.h"

// Note: For the canonical view quaternions it should be completely possible to
// determine this anaylitcally. That is the max_distance should be a
// theoretical known value
// Also: I'm not sure it matters in this case, but. We are dealing with
// quaternions on the 4d unit sphere, but measuring distance in general 4d
// space (i.e. not geodesics on the sphere). Probably something with angles
// would be better.
template <typename Q_type>
inline bool igl::snap_to_canonical_view_quat(
  const Q_type q[4],
  const Q_type threshold,
  Q_type s[4])
{
  // Copy input into output
  // CANNOT use std::copy here according to:
  // http://www.cplusplus.com/reference/algorithm/copy/
  s[0] = q[0];
  s[1] = q[1];
  s[2] = q[2];
  s[3] = q[3];

  // Normalize input quaternion
  Q_type qn[4];
  bool valid_len = 
    igl::normalize_quat(q,qn);
  // If normalizing valid then don't bother
  if(!valid_len)
  {
    return false;
  }

  // 0.290019
  const Q_type MAX_DISTANCE = 0.4;
  Q_type min_distance = 2*MAX_DISTANCE;
  int min_index = -1;
  int min_sign = 0;
  // loop over canonical view quaternions
  for(int sign = -1;sign<=1;sign+=2)
  {
    for(int i = 0; i<NUM_CANONICAL_VIEW_QUAT; i++)
    {
      Q_type distance = 0.0;
      // loop over coordinates
      for(int j = 0;j<4;j++)
      {
        distance += 
          (qn[j]-sign*igl::CANONICAL_VIEW_QUAT<Q_type>(i,j))*
          (qn[j]-sign*igl::CANONICAL_VIEW_QUAT<Q_type>(i,j));
      }
      if(min_distance > distance)
      {
        min_distance = distance;
        min_index = i;
        min_sign = sign;
      }
    }
  }

  if(MAX_DISTANCE < min_distance)
  {
    fprintf(
      stderr,
      "ERROR: found new max MIN_DISTANCE: %g\n"
      "PLEASE update snap_to_canonical_quat()",
      min_distance);
  }

  assert(min_distance < MAX_DISTANCE);
  assert(min_index >= 0);

  if( min_distance/MAX_DISTANCE <= threshold)
  {
    // loop over coordinates
    for(int j = 0;j<4;j++)
    {
      s[j] = min_sign*igl::CANONICAL_VIEW_QUAT<Q_type>(min_index,j);
    }
    return true;
  }
  return false;
}
#endif
