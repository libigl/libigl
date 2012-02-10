#include "diag.h"

#include "verbose.h"

template <typename T>
IGL_INLINE void igl::diag(
  const Eigen::SparseMatrix<T>& X, 
  Eigen::SparseVector<T>& V)
{
  // Get size of input
  int m = X.rows();
  int n = X.cols();
  V = Eigen::SparseVector<T>((m>n?n:m));

  // Iterate over outside
  for(int k=0; k<X.outerSize(); ++k)
  {
    // Iterate over inside
    for(typename Eigen::SparseMatrix<T>::InnerIterator it (X,k); it; ++it)
    {
      if(it.col() == it.row())
      {
        V.coeffRef(it.col()) += it.value();
      }
    }
  }
}

template <typename T>
IGL_INLINE void igl::diag(
  const Eigen::SparseVector<T>& V,
  Eigen::SparseMatrix<T>& X)
{
  // clear and resize output
  Eigen::DynamicSparseMatrix<T, Eigen::RowMajor> dyn_X(V.size(),V.size());

  // loop over non-zeros
  for(typename Eigen::SparseVector<T>::InnerIterator it(V); it; ++it)
  {
    dyn_X.coeffRef(it.index(),it.index()) += it.value();
  }

  X = Eigen::SparseMatrix<T>(dyn_X);
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
