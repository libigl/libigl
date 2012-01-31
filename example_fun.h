#ifndef IGL_EXAMPLE_FUN_H
#define IGL_EXAMPLE_FUN_H

#include "igl_inline.h"

namespace igl
{
  // This is an example of a function, it takes a templated parameter and
  // shovels it into cout
  //
  // Templates:
  //   T  type that supports
  // Input:
  //   input  some input of a Printable type
  // Returns true for the sake of returning something
  template <typename Printable>
  IGL_INLINE bool example_fun(const Printable & input);
}

#ifdef IGL_HEADER_ONLY
#  include "example_fun.cpp"
#endif

#endif
