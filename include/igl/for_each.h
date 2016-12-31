#ifndef IGL_FOR_EACH_H
#define IGL_FOR_EACH_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
namespace igl
{
  // FOR_EACH  Call a given function for each non-zero (i.e., explicit value
  // might actually be ==0) in a Sparse Matrix A _in order (of storage)_. This is
  // useless unless func has _side-effects_.
  //
  // Inputs:
  //   A  m by n SparseMatrix
  //   func  function handle with prototype "compatible with" `void (Index i,
  //     Index j, Scalar & v)`. Return values will be ignored.
  //
  // See also: std::for_each
  template <typename AType, typename Func>
  IGL_INLINE void for_each(
    const Eigen::SparseMatrix<AType> & A,
    const Func & func);
  template <typename DerivedA, typename Func>
  IGL_INLINE void for_each(
    const Eigen::DenseBase<DerivedA> & A,
    const Func & func);
}
#ifndef IGL_STATIC_LIBRARY
#  include "for_each.cpp"
#endif
#endif
