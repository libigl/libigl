#include "slice_into.h"

template <typename T>
IGL_INLINE void igl::slice_into(
  const Eigen::SparseMatrix<T>& X,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & C,
  Eigen::SparseMatrix<T>& Y)
{

  int xm = X.rows();
  int xn = X.cols();
  assert(R.size() == xm);
  assert(C.size() == xn);
  int ym = Y.size();
  int yn = Y.size();
  assert(R.minCoeff() >= 0);
  assert(R.maxCoeff() < ym);
  assert(C.minCoeff() >= 0);
  assert(C.maxCoeff() < yn);
  // create temporary dynamic sparse matrix
  Eigen::DynamicSparseMatrix<T, Eigen::RowMajor>  dyn_Y(Y);
  // Iterate over outside
  for(int k=0; k<X.outerSize(); ++k)
  {
    // Iterate over inside
    for(typename Eigen::SparseMatrix<T>::InnerIterator it (X,k); it; ++it)
    {
      dyn_Y.coeffRef(R(it.row()),C(it.col())) = it.value();
    }
  }
  Y = Eigen::SparseMatrix<T>(dyn_Y);
}

template <typename T, const int W, const int H>
IGL_INLINE void igl::slice_into(
  const Eigen::Matrix<T,W,H> & X,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & C,
  Eigen::Matrix<T,W,H> & Y)
{

  int xm = X.rows();
  int xn = X.cols();
  assert(R.size() == xm);
  assert(C.size() == xn);
  int ym = Y.size();
  int yn = Y.size();
  assert(R.minCoeff() >= 0);
  assert(R.maxCoeff() < ym);
  assert(C.minCoeff() >= 0);
  assert(C.maxCoeff() < yn);

  // Build reindexing maps for columns and rows, -1 means not in map
  Eigen::Matrix<int,Eigen::Dynamic,1> RI;
  RI.resize(xm);
  for(int i = 0;i<xm;i++)
  {
    for(int j = 0;j<xn;j++)
    {
      Y(R(i),C(j)) = X(i,j);
    }
  }
}

template <typename T>
IGL_INLINE void igl::slice_into(
  const Eigen::Matrix<T,Eigen::Dynamic,1> & X,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
  Eigen::Matrix<T,Eigen::Dynamic,1> & Y)
{
  // phony column indices
  Eigen::Matrix<int,Eigen::Dynamic,1> C;
  C.resize(1);
  C(0) = 0;
  return igl::slice_into(X,R,C,Y);
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
