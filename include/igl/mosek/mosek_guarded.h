#ifndef IGL_MOSEK_GUARDED_H
#define IGL_MOSEK_GUARDED_H
#include "../igl_inline.h"

#include "mosek.h"
namespace igl
{
  // Little function to wrap around mosek call to handle errors
  // 
  // Inputs:
  //   r  mosek error code returned from mosek call
  // Returns r untouched
  IGL_INLINE MSKrescodee mosek_guarded(const MSKrescodee r);
}

#ifdef IGL_HEADER_ONLY
#  include "mosek_guarded.cpp"
#endif

#endif

