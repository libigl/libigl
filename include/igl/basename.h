#ifndef IGL_BASENAME_H
#define IGL_BASENAME_H
#include "igl_inline.h"

#include <string>

namespace igl
{
  // Function like PHP's basename
  // Input:
  //  path  string containing input path
  // Returns string containing basename (see php's basename)
  //
  // See also: dirname, pathinfo
  IGL_INLINE std::string basename(const std::string & path);
}

#ifdef IGL_HEADER_ONLY
#  include "basename.cpp"
#endif

#endif
