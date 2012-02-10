#ifndef IGL_FILE_EXISTS_H
#define IGL_FILE_EXISTS_H
#include "igl_inline.h"
namespace igl
{
  // Check if a file or directory exists like PHP's file_exists function:
  // http://php.net/manual/en/function.file-exists.php
  // Input:
  //   filename  path to file
  // Returns true if file exists and is readable and false if file doesn't
  // exist or *is not readable*
  IGL_INLINE bool file_exists(const char * filename);
}

#ifdef IGL_HEADER_ONLY
#  include "file_exists.cpp"
#endif

#endif
