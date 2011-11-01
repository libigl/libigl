#ifndef IGL_MAX_SIZE_H
#define IGL_MAX_SIZE_H
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
  inline int max_size(const std::vector<T> & V);
}

// implementation 

template <typename T>
inline int igl::max_size(const std::vector<T> & V)
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
#endif
