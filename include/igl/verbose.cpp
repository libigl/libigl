#include "verbose.h"

// http://channel9.msdn.com/forums/techoff/254707-wrapping-printf-in-c/
IGL_INLINE int igl::verbose(const char * msg,...)
{
#ifdef VERBOSE
  va_list argList;
  va_start(argList, msg);
  int count = vprintf(msg, argList);
  va_end(argList);
  return count;
#else
  return 0;
#endif
}
