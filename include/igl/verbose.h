#ifndef IGL_VERBOSE_H
#define IGL_VERBOSE_H
#include "igl_inline.h"
// Provide a macro for printf, called verbose that functions exactly like
// printf if VERBOSE is defined and does exactly nothing if VERBOSE is
// undefined
#include <cstdio>
#ifdef VERBOSE
#  include <cstdarg>
#endif

namespace igl
{
  IGL_INLINE int verbose(const char * msg,...);
}

#ifdef IGL_HEADER_ONLY
#  include "verbose.cpp"
#endif

#endif
