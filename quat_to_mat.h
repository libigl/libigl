#ifndef IGL_QUAT_TO_MAT_H
#define IGL_QUAT_TO_MAT_H
#include "igl_inline.h"
// Name history:
//   quat2mat  until 16 Sept 2011
namespace igl
{
  // Convert a quaternion to a 4x4 matrix
  // A Quaternion, q, is defined here as an arrays of four scalars (x,y,z,w),
  // such that q = x*i + y*j + z*k + w
  // Input:
  //   quat  pointer to four elements of quaternion (x,y,z,w)  
  // Output:
  //   mat  pointer to 16 elements of matrix
  template <typename Q_type>
  IGL_INLINE void quat_to_mat(const Q_type * quat, Q_type * mat);
}

#ifdef IGL_HEADER_ONLY
#  include "quat_to_mat.cpp"
#endif

#endif
