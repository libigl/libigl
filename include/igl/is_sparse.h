#ifndef IGL_IS_SPARSE_H
#define IGL_IS_SPARSE_H
#include "igl_inline.h"
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>
namespace igl
{
  // Determine if a matrix A is sparse
  //
  // Template:
  //   T,DerivedA defines scalar type
  // Inputs:
  //   A  matrix in question
  // Returns true if A is represented with a sparse matrix
  template <typename T>
  IGL_INLINE bool is_sparse(
    const Eigen::SparseMatrix<T> & A);
  template <typename DerivedA>
  IGL_INLINE bool is_sparse(
    const Eigen::PlainObjectBase<DerivedA>& A);
}

#ifdef IGL_HEADER_ONLY
#  include "is_sparse.cpp"
#endif

#endif

