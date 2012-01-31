#include "repmat.h"

template <typename T,const int W, const int H>
IGL_INLINE void igl::repmat(
  const Eigen::Matrix<T,W,H> & A,
  const int r,
  const int c,
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & B)
{
  assert(r>0);
  assert(c>0);
  // Make room for output
  B.resize(r*A.rows(),c*A.cols());

  // copy tiled blocks
  for(int i = 0;i<r;i++)
  {
    for(int j = 0;j<c;j++)
    {
      B.block(i*A.rows(),j*A.cols(),A.rows(),A.cols()) = A;
    }
  }
}

template <typename T>
IGL_INLINE void igl::repmat(
  const Eigen::SparseMatrix<T> & A,
  const int r,
  const int c,
  Eigen::SparseMatrix<T> & B)
{
  assert(r>0);
  assert(c>0);
  B.resize(r*A.rows(),c*A.cols());
  B.reserve(r*c*A.nonZeros());
  for(int i = 0;i<r;i++)
  {
    for(int j = 0;j<c;j++)
    {
      // Loop outer level
      for (int k=0; k<A.outerSize(); ++k)
      {
        // loop inner level
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(A,k); it; ++it)
        {
          B.insert(i*A.rows()+it.row(),j*A.cols() + it.col()) = it.value();
        }
      }
    }
  }
  B.finalize();
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
