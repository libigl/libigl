#ifndef IGL_REPDIAG_H
#define IGL_REPDIAG_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace igl
{
  // REPDIAG repeat a matrix along the diagonal a certain number of times, so
  // that if A is a m by n matrix and we want to repeat along the diagonal d
  // times, we get a m*d by n*d matrix B such that:
  // B( (k*m+1):(k*m+1+m-1), (k*n+1):(k*n+1+n-1)) = A 
  // for k from 0 to d-1
  //
  // Inputs:
  //   A  m by n matrix we are repeating along the diagonal. May be dense or
  //     sparse
  //   d  number of times to repeat A along the diagonal
  // Outputs:
  //   B  m*d by n*d matrix with A repeated d times along the diagonal,
  //     will be dense or sparse to match A
  //

  // Sparse version
  template <typename T>
  inline void repdiag(
    const Eigen::SparseMatrix<T>& A,
    const int d,
    Eigen::SparseMatrix<T>& B);
  // Dense version
  template <typename T>
  inline void repdiag(
    const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & A,
    const int d,
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & B);
}

// Implementation

template <typename T>
inline void igl::repdiag(
  const Eigen::SparseMatrix<T>& A,
  const int d,
  Eigen::SparseMatrix<T>& B)
{
  int m = A.rows();
  int n = A.cols();

  Eigen::DynamicSparseMatrix<T, Eigen::RowMajor> 
    dyn_B(m*d,n*d);

  // loop over reps
  for(int i=0;i<d;i++)
  {
    // Loop outer level
    for (int k=0; k<A.outerSize(); ++k)
    {
      // loop inner level
      for (typename Eigen::SparseMatrix<T>::InnerIterator it(A,k); it; ++it)
      {
        dyn_B.coeffRef(i*m+it.row(),i*n+it.col()) += it.value();
      }
    }
  }

  B = Eigen::SparseMatrix<T>(dyn_B);
}

template <typename T>
inline void igl::repdiag(
  const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & A,
  const int d,
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & B)
{
  int m = A.rows();
  int n = A.cols();
  B.resize(m*d,n*d);
  B.array() *= 0;
  for(int i = 0;i<d;i++)
  {
    B.block(i*m,i*n,m,n) = A;
  }
}

#endif
