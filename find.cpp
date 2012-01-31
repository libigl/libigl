#include "find.h"

#include "verbose.h"
  
template <typename T>
IGL_INLINE void igl::find(
  const Eigen::SparseMatrix<T>& X,
  Eigen::Matrix<int,Eigen::Dynamic,1> & I,
  Eigen::Matrix<int,Eigen::Dynamic,1> & J,
  Eigen::Matrix<T,Eigen::Dynamic,1> & V)
{
  // Resize outputs to fit nonzeros
  I.resize(X.nonZeros());
  J.resize(X.nonZeros());
  V.resize(X.nonZeros());

  int i = 0;
  // Iterate over outside
  for(int k=0; k<X.outerSize(); ++k)
  {
    // Iterate over inside
    for(typename Eigen::SparseMatrix<T>::InnerIterator it (X,k); it; ++it)
    {
      V(i) = it.value();
      I(i) = it.row();
      J(i) = it.col();
      i++;
    }
  }
}
  
template <typename T>
IGL_INLINE void igl::find(
  const Eigen::SparseVector<T>& X,
  Eigen::Matrix<int,Eigen::Dynamic,1> & I,
  Eigen::Matrix<T,Eigen::Dynamic,1> & V)
{
  // Resize outputs to fit nonzeros
  I.resize(X.nonZeros());
  V.resize(X.nonZeros());

  int i = 0;
  // loop over non-zeros
  for(typename Eigen::SparseVector<T>::InnerIterator it(X); it; ++it)
  {
    I(i) = it.index();
    V(i) = it.value();
    i++;
  }
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
