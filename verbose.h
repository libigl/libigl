#ifndef IGL_VERBOSE_H
#define IGL_VERBOSE_H
// Provide a macro for printf, called verbose that functions exactly like
// printf if VERBOSE is defined and does exactly nothing if VERBOSE is
// undefined
#include <cstdio>
#ifdef VERBOSE
#  include <cstdarg>
#endif

namespace igl
{
  inline int verbose(const char * msg,...);
};

// Implementation

// http://channel9.msdn.com/forums/techoff/254707-wrapping-printf-in-c/
inline int igl::verbose(const char * msg,...)
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
#endif
