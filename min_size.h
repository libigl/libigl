#ifndef IGL_MIN_SIZE_H
#define IGL_MIN_SIZE_H
#include <vector>

namespace igl
{
  // Determine min size of lists in a vector
  // Template:
  //   T  some list type object that implements .size()
  // Inputs:
  //   V  vector of list types T
  // Returns min .size() found in V, returns -1 if V is empty
  template <typename T>
  inline int min_size(const std::vector<T> & V);
}

// implementation 

template <typename T>
inline int igl::min_size(const std::vector<T> & V)
{
  int min_size = -1;
  for(
    typename std::vector<T>::const_iterator iter = V.begin();
    iter != V.end(); 
    iter++)
  {
    int size = (int)iter->size();
    // have to handle base case
    if(min_size == -1)
    {
      min_size = size;
    }else{
      min_size = (min_size < size ? min_size : size);
    }
  }
  return min_size;
}
#endif
