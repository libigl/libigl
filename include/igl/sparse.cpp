#include "sparse.h"

template <class IndexVector, class ValueVector, typename T>
IGL_INLINE void igl::sparse(
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
IGL_INLINE void igl::sparse(
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

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
