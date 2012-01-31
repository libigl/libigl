#ifndef IGL_IS_WRITABLE_H
#define IGL_IS_WRITABLE_H
#include "igl_inline.h"
namespace igl
{
  // Check if a file exists *and* is writable like PHP's is_writable function:
  // http://www.php.net/manual/en/function.is-writable.php
  // Input:
  //   filename  path to file
  // Returns true if file exists and is writable and false if file doesn't
  // exist or *is not writable*
  //
  // Note: Windows version will not test group and user id
  IGL_INLINE bool is_writable(const char * filename);
}

#ifdef IGL_HEADER_ONLY
#  include "is_writable.cpp"
#endif

#endif
