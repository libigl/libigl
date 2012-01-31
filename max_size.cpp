#include "max_size.h"


template <typename T>
IGL_INLINE int igl::max_size(const std::vector<T> & V)
{
  int max_size = -1;
  for(
    typename std::vector<T>::const_iterator iter = V.begin();
    iter != V.end(); 
    iter++)
  {
    int size = (int)iter->size();
    max_size = (max_size > size ? max_size : size);
  }
  return max_size;
}


#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
