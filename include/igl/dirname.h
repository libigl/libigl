#ifndef IGL_DIRNAME_H
#define IGL_DIRNAME_H
#include "igl_inline.h"

#include <string>

namespace igl
{
  // Function like PHP's dirname
  // Input:
  //  path  string containing input path
  // Returns string containing dirname (see php's dirname)
  //
  // See also: basename, pathinfo
  IGL_INLINE std::string dirname(const std::string & path);
}

#ifdef IGL_HEADER_ONLY
#  include "dirname.cpp"
#endif

#endif
