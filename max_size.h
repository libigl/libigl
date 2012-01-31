#ifndef IGL_MAX_SIZE_H
#define IGL_MAX_SIZE_H
#include "igl_inline.h"
#include <vector>

namespace igl
{
  // Determine max size of lists in a vector
  // Template:
  //   T  some list type object that implements .size()
  // Inputs:
  //   V  vector of list types T
  // Returns max .size() found in V, returns -1 if V is empty
  template <typename T>
  IGL_INLINE int max_size(const std::vector<T> & V);
}

#ifdef IGL_HEADER_ONLY
#  include "max_size.cpp"
#endif

#endif
