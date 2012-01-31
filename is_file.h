#ifndef IGL_IS_FILE_H
#define IGL_IS_FILE_H
#include "igl_inline.h"
namespace igl
{
  // Act like php's is_file function
  // http://php.net/manual/en/function.is-file.php
  // Tells whether the given filename is a regular file.
  // Input:
  //   filename  Path to the file. If filename is a relative filename, it will
  //     be checked relative to the current working directory. 
  // Returns TRUE if the filename exists and is a regular file, FALSE
  // otherwise.
  IGL_INLINE bool is_file(const char * filename);

}

#ifdef IGL_HEADER_ONLY
#  include "is_file.cpp"
#endif

#endif
