#include "invert_diag.h"
#include "diag.h"

template <typename T>
IGL_INLINE void igl::invert_diag(
  const Eigen::SparseMatrix<T>& X, 
  Eigen::SparseMatrix<T>& Y)
{
#ifndef NDEBUG
  typename Eigen::SparseVector<T> dX;
  igl::diag(X,dX);
  // Check that there are no zeros along the diagonal
  assert(dX.nonZeros() == dX.size());
#endif
  // http://www.alecjacobson.com/weblog/?p=2552
  if(&Y != &X)
  {
    Y = X;
  }
  // Iterate over outside
  for(int k=0; k<Y.outerSize(); ++k)
  {
    // Iterate over inside
    for(typename Eigen::SparseMatrix<T>::InnerIterator it (Y,k); it; ++it)
    {
      if(it.col() == it.row())
      {
        T v = it.value();
        assert(v != 0);
        v = ((T)1.0)/v;
        Y.coeffRef(it.row(),it.col()) = v;
      }
    }
  }
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template void igl::invert_diag<double>(Eigen::SparseMatrix<double, 0, int> const&, Eigen::SparseMatrix<double, 0, int>&);
#endif
