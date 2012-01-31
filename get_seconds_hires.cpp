#include "get_seconds_hires.h"

#if _WIN32
#  include <windows.h>
IGL_INLINE double igl::get_seconds_hires()
{
  LARGE_INTEGER li_freq, li_current;
  const bool ret = QueryPerformanceFrequency(&li_freq);
  const bool ret2 = QueryPerformanceCounter(&li_current);
  assert(ret && ret2);
  assert(li_freq.QuadPart > 0);
  return double(li_current.QuadPart) / double(li_freq.QuadPart);
}
#else
#  include "get_seconds.h"
IGL_INLINE double igl::get_seconds_hires()
{
  // Sorry I've no idea how performance counters work on Mac...
  return igl::get_seconds();
}
#endif
