#ifndef IGL_REDUX_H
#define IGL_REDUX_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
namespace igl
{
  // REDUX Perform reductions on the rows or columns of a SparseMatrix. This is
  // _similar_ to DenseBase::redux, but different in two important ways:
  //  1. (unstored) Zeros are **not** "visited", however if the first element
  //     in the column/row  does not appear in the first row/column then the
  //     reduction is assumed to start with zero. In this way, "any", "all",
  //     "count"(non-zeros) work as expected. This means it is **not** possible
  //     to use this to count (implicit) zeros.
  //  2. This redux is more powerful in the sense that A and B may have
  //     different types. This makes it possible to count the number of
  //     non-zeros in a SparseMatrix<bool> A into a VectorXi B.
  //
  // Inputs:
  //   A  m by n sparse matrix
  //   dim  dimension along which to sum (1 or 2)
  //   func  function handle with the prototype `X(Y a, I i, J j, Z b)` where a
  //     is the running value, b is A(i,j)
  // Output:
  //   S  n-long sparse vector (if dim == 1) 
  //   or
  //   S  m-long sparse vector (if dim == 2)
  template <typename AType, typename Func, typename DerivedB>
  IGL_INLINE void redux(
    const Eigen::SparseMatrix<AType> & A,
    const int dim,
    const Func & func,
    Eigen::PlainObjectBase<DerivedB> & B);
}
#ifndef IGL_STATIC_LIBRARY
#  include "redux.cpp"
#endif
#endif
