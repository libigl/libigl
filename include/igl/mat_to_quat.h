#ifndef IGL_MAT_TO_QUAT_H
#define IGL_MAT_TO_QUAT_H
#include "igl_inline.h"
namespace igl
{
  // Convert a OpenGL (rotation) matrix to a quaternion
  //
  // Input:
  //   m  16-element opengl rotation matrix
  // Output:
  //   q  4-element  quaternion (not normalized)
  template <typename Q_type>
  IGL_INLINE void mat4_to_quat(const Q_type * m, Q_type * q);
  // Input:
  //   m  9-element opengl rotation matrix
  template <typename Q_type>
  IGL_INLINE void mat3_to_quat(const Q_type * m, Q_type * q);
}

#ifdef IGL_HEADER_ONLY
#  include "mat_to_quat.cpp"
#endif

#endif

