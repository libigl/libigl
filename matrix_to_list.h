#ifndef IGL_MATRIX_TO_LIST_H
#define IGL_MATRIX_TO_LIST_H
#include "igl_inline.h"
#include <vector>
namespace igl
{
  // Convert a matrix to a list (std::vector) of row vectors of the same size
  // Template: 
  //   Mat  Matrix type, must implement:
  //     .resize(m,n)
  //     .row(i) = Row
  //   T  type that can be safely cast to type in Mat via '='
  // Inputs:
  //   V  a m-long list of vectors of size n
  // Outputs:
  //   M  an m by n matrix
  //
  // See also: list_to_matrix
  template <typename Mat, typename T>
  IGL_INLINE void matrix_to_list(
    const Mat & M, 
    std::vector<std::vector<T > > & V);
}

#ifdef IGL_HEADER_ONLY
#  include "matrix_to_list.cpp"
#endif

#endif

