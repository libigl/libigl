#include "dot.h"

// http://www.antisphere.com/Wiki/tools:anttweakbar
IGL_INLINE double igl::dot(
  const double *a, 
  const double *b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
