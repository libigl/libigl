#ifndef IGL_NORMALIZE_ROW_SUMS_H
#define IGL_NORMALIZE_ROW_SUMS_H
#include "igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  // Normalize the rows in A so that their sums are each 1 and place the new
  // entries in B
  // Inputs:
  //   A  #rows by k input matrix
  // Outputs:
  //   B  #rows by k input matrix, can be the same as A
  template <typename DerivedA, typename DerivedB>
  IGL_INLINE void normalize_row_sums(
    const Eigen::MatrixBase<DerivedA>& A,
    Eigen::MatrixBase<DerivedB> & B);
}

#ifdef IGL_HEADER_ONLY
#  include "normalize_row_sums.cpp"
#endif

#endif
