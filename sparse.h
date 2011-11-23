#ifndef IGL_SPARSE_H
#define IGL_SPARSE_H
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>
namespace igl
{
  // Build a sparse matrix from list of indices and values (I,J,V), functions
  // like the sparse function in matlab
  //
  // Templates:
  //   IndexVector  list of indices, value should be non-negative and should
  //     expect to be cast to an index. Must implement operator(i) to retrieve
  //     ith element
  //   ValueVector  list of values, value should be expect to be cast to type
  //     T. Must implement operator(i) to retrieve ith element
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Input:
  //   I  nnz vector of row indices of non zeros entries in X
  //   J  nnz vector of column indices of non zeros entries in X
  //   V  nnz vector of non-zeros entries in X
  //   Optional:
  //     m  number of rows
  //     n  number of cols
  // Outputs:
  //   X  m by n matrix of type T whose entries are to be found 
  //
  template <class IndexVector, class ValueVector, typename T>
  inline void sparse(
    const IndexVector & I,
    const IndexVector & J,
    const ValueVector & V,
    Eigen::SparseMatrix<T>& X);
  template <class IndexVector, class ValueVector, typename T>
  inline void sparse(
    const IndexVector & I,
    const IndexVector & J,
    const ValueVector & V,
    const size_t m,
    const size_t n,
    Eigen::SparseMatrix<T>& X);
}

// Implementation

template <class IndexVector, class ValueVector, typename T>
inline void igl::sparse(
  const IndexVector & I,
  const IndexVector & J,
  const ValueVector & V,
  Eigen::SparseMatrix<T>& X)
{
  size_t m = (size_t)I.maxCoeff()+1;
  size_t n = (size_t)J.maxCoeff()+1;
  return igl::sparse(I,J,V,m,n,X);
}

#include "verbose.h"
template <class IndexVector, class ValueVector, typename T>
inline void igl::sparse(
  const IndexVector & I,
  const IndexVector & J,
  const ValueVector & V,
  const size_t m,
  const size_t n,
  Eigen::SparseMatrix<T>& X)
{
  assert((int)I.maxCoeff() < (int)m);
  assert((int)I.minCoeff() >= 0);
  assert((int)J.maxCoeff() < (int)n);
  assert((int)J.minCoeff() >= 0);
  assert(I.size() == J.size());
  assert(J.size() == V.size());
  // Really we just need .size() to be the same, but this is safer
  assert(I.rows() == J.rows());
  assert(J.rows() == V.rows());
  assert(I.cols() == J.cols());
  assert(J.cols() == V.cols());
  // number of values
  int nv = V.size();

  Eigen::DynamicSparseMatrix<T, Eigen::RowMajor> dyn_X(m,n);
  // over estimate the number of entries
  dyn_X.reserve(I.size());
  for(int i = 0;i < nv;i++)
  {
    dyn_X.coeffRef((int)I(i),(int)J(i)) += (T)V(i);
  }
  X = Eigen::SparseMatrix<T>(dyn_X);
}

#endif
