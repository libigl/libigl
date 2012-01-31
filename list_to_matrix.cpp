#include "list_to_matrix.h"

#include <cassert>
#include <cstdio>

#include "max_size.h"
#include "min_size.h"
#define VERBOSE
#include "verbose.h"

template <typename T, class Mat>
IGL_INLINE bool igl::list_to_matrix(const std::vector<std::vector<T > > & V,Mat & M)
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

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
