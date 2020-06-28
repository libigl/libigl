#ifndef IGL_ALWAYS_TRY_NEVER_CARE_H
#define IGL_ALWAYS_TRY_NEVER_CARE_H
#include "igl_inline.h"
#include "decimate_func_types.h"
namespace igl
{
  // Outputs:
  //   always_try  function that always returns true
  //   never_care  fuction that is always a no-op
  IGL_INLINE void always_try_never_care(
    decimate_pre_collapse_func  & always_try,
    decimate_post_collapse_func & never_care);
};

#ifndef IGL_STATIC_LIBRARY
#  include "always_try_never_care.cpp"
#endif

#endif 

