#ifndef IGL_GET_SECONDS_H
#define IGL_GET_SECONDS_H
#include "igl_inline.h"

namespace igl
{
  // Return the current time in seconds since program start
  IGL_INLINE double get_seconds();

}

#ifdef IGL_HEADER_ONLY
#  include "get_seconds.cpp"
#endif

#endif
