#ifndef IGL_SLICE_H
#define IGL_SLICE_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>

namespace igl
{
  // Act like the matlab X(row_indices,col_indices) operator
  // 
  // Inputs:
  //   X  m by n matrix
  //   R  list of row indices
  //   C  list of column indices
  // Output:
  //   Y  #R by #C matrix
  template <typename T>
  inline void slice(
    const Eigen::SparseMatrix<T>& X,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & C,
    Eigen::SparseMatrix<T>& Y);

  template <typename T, const int W, const int H>
  inline void slice(
    const Eigen::Matrix<T,W,H> & X,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & C,
    Eigen::Matrix<T,W,H> & Y);

  template <typename T>
  inline void slice(
    const Eigen::Matrix<T,Eigen::Dynamic,1> & X,
    const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
    Eigen::Matrix<T,Eigen::Dynamic,1> & Y);
}

// Implementation

template <typename T>
inline void igl::slice(
  const Eigen::SparseMatrix<T>& X,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & C,
  Eigen::SparseMatrix<T>& Y)
{
  int xm = X.rows();
  int xn = X.cols();
  int ym = R.size();
  int yn = C.size();
  assert(R.minCoeff() >= 0);
  assert(R.maxCoeff() < xm);
  assert(C.minCoeff() >= 0);
  assert(C.maxCoeff() < xn);
  // Build reindexing maps for columns and rows, -1 means not in map
  Eigen::Matrix<int,Eigen::Dynamic,1> RI;
  RI.resize(xm);
  // initialize to -1
  for(int i = 0;i<xm;i++)
  {
    RI(i) = -1;
  }
  for(int i = 0;i<ym;i++)
  {
    RI(R(i)) = i;
  }
  Eigen::Matrix<int,Eigen::Dynamic,1> CI;
  CI.resize(xn);
  // initialize to -1
  for(int i = 0;i<xn;i++)
  {
    CI(i) = -1;
  }
  for(int i = 0;i<yn;i++)
  {
    CI(C(i)) = i;
  }
  // Resize output
  Eigen::DynamicSparseMatrix<T, Eigen::RowMajor> 
    dyn_Y(ym,yn);
  // Take a guess at the number of nonzeros (this assumes uniform distribution
  // not banded or heavily diagonal)
  dyn_Y.reserve((X.nonZeros()/(X.rows()*X.cols())) * (ym*yn));
  // Iterate over outside
  for(int k=0; k<X.outerSize(); ++k)
  {
    // Iterate over inside
    for(typename Eigen::SparseMatrix<T>::InnerIterator it (X,k); it; ++it)
    {
      if(RI(it.row()) >= 0 && CI(it.col()) >= 0)
      {
        dyn_Y.coeffRef(RI(it.row()),CI(it.col())) = it.value();
      }
    }
  }
  Y = Eigen::SparseMatrix<T>(dyn_Y);
}

template <typename T, const int W, const int H>
inline void igl::slice(
  const Eigen::Matrix<T,W,H> & X,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & C,
  Eigen::Matrix<T,W,H> & Y)
{
  int xm = X.rows();
  int xn = X.cols();
  int ym = R.size();
  int yn = C.size();
  assert(R.minCoeff() >= 0);
  assert(R.maxCoeff() < xm);
  assert(C.minCoeff() >= 0);
  assert(C.maxCoeff() < xn);
  // Build reindexing maps for columns and rows, -1 means not in map
  Eigen::Matrix<int,Eigen::Dynamic,1> RI;
  RI.resize(xm);
  // initialize to -1
  for(int i = 0;i<xm;i++)
  {
    RI(i) = -1;
  }
  for(int i = 0;i<ym;i++)
  {
    RI(R(i)) = i;
  }
  Eigen::Matrix<int,Eigen::Dynamic,1> CI;
  CI.resize(xn);
  // initialize to -1
  for(int i = 0;i<xn;i++)
  {
    CI(i) = -1;
  }
  for(int i = 0;i<yn;i++)
  {
    CI(C(i)) = i;
  }
  // Resize output
  Y.resize(ym,yn);
  for(int i = 0;i<xm;i++)
  {
    for(int j = 0;j<xn;j++)
    {
      if(RI(i) >= 0 && CI(j) >= 0)
      {
        Y(RI(i),CI(j)) = X(i,j);
      }
    }
  }
}

template <typename T>
inline void igl::slice(
  const Eigen::Matrix<T,Eigen::Dynamic,1> & X,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
  Eigen::Matrix<T,Eigen::Dynamic,1> & Y)
{
  // phony column indices
  Eigen::Matrix<int,Eigen::Dynamic,1> C;
  C.resize(1);
  C(0) = 0;
  return igl::slice(X,R,C,Y);
}

#endif

