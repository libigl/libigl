#ifndef IGL_LIST_TO_MATRIX_H
#define IGL_LIST_TO_MATRIX_H
#include <vector>
namespace igl
{
  // Convert a list (std::vector) of row vectors of the same length to a matrix
  // Template: 
  //   T  type that can be safely cast to type in Mat via '='
  //   Mat  Matrix type, must implement:
  //     .resize(m,n)
  //     .row(i) = Row
  // Inputs:
  //   V  a m-long list of vectors of size n
  // Outputs:
  //   M  an m by n matrix
  // Returns true on success, false on errors
  template <typename T, class Mat>
  inline bool list_to_matrix(const std::vector<std::vector<T > > & V,Mat & M);
}

// Implementation
#include <cassert>
#include <cstdio>

#include "max_size.h"
#include "min_size.h"
#define VERBOSE
#include "verbose.h"

template <typename T, class Mat>
inline bool igl::list_to_matrix(const std::vector<std::vector<T > > & V,Mat & M)
{
  // number of columns
  int m = V.size();
  if(m == 0)
  {
    fprintf(stderr,"Error: list_to_matrix() list is empty()\n");
    return false;
  }
  // number of rows
  int n = igl::min_size(V);
  if(n != igl::max_size(V))
  {
    fprintf(stderr,"Error: list_to_matrix()"
      " list elements are not all the same size\n");
    return false;
  }
  assert(n != -1);
  // Resize output
  M.resize(m,n);

  // Loop over rows
  for(int i = 0;i<m;i++)
  {
    // Loop over cols
    for(int j = 0;j<n;j++)
    {
      M(i,j) = V[i][j];
    }
  }

  return true;
}
#endif
