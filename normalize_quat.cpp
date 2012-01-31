#include "normalize_quat.h"

#include "EPS.h"

template <typename Q_type>
IGL_INLINE bool igl::normalize_quat(
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

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
