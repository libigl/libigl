#ifndef IGL_NORMALIZE_ROWS_H
#define IGL_NORMALIZE_ROWS_H
#include "igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  // Normalize the rows in A so that their lengths are each 1 and place the new
  // entries in B
  // Inputs:
  //   A  #rows by k input matrix
  // Outputs:
  //   B  #rows by k input matrix, can be the same as A
  IGL_INLINE void normalize_rows(const Eigen::MatrixXd & A, Eigen::MatrixXd & B);
}

#ifdef IGL_HEADER_ONLY
#  include "normalize_rows.cpp"
#endif

#endif
