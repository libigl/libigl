#ifndef IGL_ROWS_TO_MATRIX_H
#define IGL_ROWS_TO_MATRIX_H
#include <vector>
namespace igl
{
  // Convert a list (std::vector) of row vectors of the same length to a matrix
  // Template: 
  //   Row  row vector type, must implement:
  //     .size()
  //   Mat  Matrix type, must implement:
  //     .resize(m,n)
  //     .row(i) = Row
  // Inputs:
  //   V  a m-long list of vectors of size n
  // Outputs:
  //   M  an m by n matrix
  // Returns true on success, false on errors
  template <class Row, class Mat>
  bool rows_to_matrix(const std::vector<Row> & V,Mat & M);
}

// Implementation
#include <cassert>
#include <cstdio>

#include "max_size.h"
#include "min_size.h"

template <class Row, class Mat>
bool igl::rows_to_matrix(const std::vector<Row> & V,Mat & M)
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
  int i = 0;
  typename std::vector<Row>::const_iterator iter = V.begin();
  while(iter != V.end())
  {
    M.row(i) = V[i];
    // increment index and iterator
    i++;
    iter++;
  }

  return true;
}
#endif
