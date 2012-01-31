#ifndef IGL_IS_READABLE_H
#define IGL_IS_READABLE_H
#include "igl_inline.h"
namespace igl
{
  // Check if a file is reabable like PHP's is_readable function:
  // http://www.php.net/manual/en/function.is-readable.php
  // Input:
  //   filename  path to file
  // Returns true if file exists and is readable and false if file doesn't
  // exist or *is not readable*
  //
  // Note: Windows version will not check user or group ids
  IGL_INLINE bool is_readable(const char * filename);
}

#ifdef IGL_HEADER_ONLY
#  include "is_readable.cpp"
#endif

#endif
