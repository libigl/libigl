#ifndef IGL_GET_SECONDS_HIRES_H
#define IGL_GET_SECONDS_HIRES_H
#include "igl_inline.h"

namespace igl
{
  // Return the current time in seconds using performance counters
  IGL_INLINE double get_seconds_hires();
}

#ifdef IGL_HEADER_ONLY
#  include "get_seconds_hires.cpp"
#endif

#endif
