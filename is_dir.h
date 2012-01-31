#ifndef IGL_IS_DIR_H
#define IGL_IS_DIR_H
#include "igl_inline.h"
namespace igl
{
  // Act like php's is_dir function
  // http://php.net/manual/en/function.is-dir.php
  // Tells whether the given filename is a directory.
  // Input:
  //   filename  Path to the file. If filename is a relative filename, it will
  //     be checked relative to the current working directory. 
  // Returns TRUE if the filename exists and is a directory, FALSE
  // otherwise.
  IGL_INLINE bool is_dir(const char * filename);

}

#ifdef IGL_HEADER_ONLY
#  include "is_dir.cpp"
#endif

#endif
