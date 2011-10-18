#ifndef IGL_REPMAT_H
#define IGL_REPMAT_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace igl
{
  // Ideally this is a super overloaded function that behaves just like
  // matlab's repmat

  // Replicate and tile a matrix
  //
  // Templates:
  //   T  should be a eigen matrix primitive type like int or double
  // Inputs:
  //   A  m by n input matrix
  //   r  number of row-direction copies
  //   c  number of col-direction copies
  // Outputs:
  //   B  r*m by c*n output matrix
  //
  template <typename T,const int W, const int H>
  inline void repmat(
    const Eigen::Matrix<T,W,H> & A,
    const int r,
    const int c,
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & B);
  template <typename T>
  inline void repmat(
    const Eigen::SparseMatrix<T> & A,
    const int r,
    const int c,
    Eigen::SparseMatrix<T> & B);
}

// Implementation

template <typename T,const int W, const int H>
inline void igl::repmat(
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
inline void igl::repmat(
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
#endif
