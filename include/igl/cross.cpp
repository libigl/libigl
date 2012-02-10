#include "cross.h"

// http://www.antisphere.com/Wiki/tools:anttweakbar
IGL_INLINE void igl::cross(
  const double *a, 
  const double *b,
  double *out)
{
  out[0] = a[1]*b[2]-a[2]*b[1];
  out[1] = a[2]*b[0]-a[0]*b[2];
  out[2] = a[0]*b[1]-a[1]*b[0];
}
