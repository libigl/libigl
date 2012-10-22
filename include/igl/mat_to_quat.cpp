#include "mat_to_quat.h"
#include <cmath>

// This could be replaced by something fast
template <typename Q_type>
static inline Q_type ReciprocalSqrt( const Q_type x )
{
  return 1.0/sqrt(x);
}

// Converts row major order matrix to quat
// http://software.intel.com/sites/default/files/m/d/4/1/d/8/293748.pdf
template <typename Q_type>
IGL_INLINE void igl::mat4_to_quat(const Q_type * m, Q_type * q)
{
  Q_type t = + m[0 * 4 + 0] + m[1 * 4 + 1] + m[2 * 4 + 2] + 1.0f; 
  Q_type s = ReciprocalSqrt( t ) * 0.5f;
  q[3] = s * t;
  q[2] = ( m[0 * 4 + 1] - m[1 * 4 + 0] ) * s; 
  q[1] = ( m[2 * 4 + 0] - m[0 * 4 + 2] ) * s; 
  q[0] = ( m[1 * 4 + 2] - m[2 * 4 + 1] ) * s;
}


#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template void igl::mat4_to_quat<double>(double const*, double*);
#endif
