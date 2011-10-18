#ifndef IGL_TRANSPOSE_BLOCKS
#define IGL_TRANSPOSE_BLOCKS

#include <Eigen/Core>

namespace igl
{
  // Templates:
  //   T  should be a eigen matrix primitive type like int or double
  // Inputs:
  //   A  m*k by n (dim: 1) or m by n*k (dim: 2) eigen Matrix of type T values
  //   k  number of blocks
  //   dim  dimension in which to transpose
  // Output
  //   B  n*k by m (dim: 1) or n by m*k (dim: 2) eigen Matrix of type T values,
  //   NOT allowed to be the same as A
  //
  // Example:
  // A = [
  //   1   2   3   4
  //   5   6   7   8
  // 101 102 103 104
  // 105 106 107 108
  // 201 202 203 204
  // 205 206 207 208];
  // transpose_blocks(A,1,3,B);
  // B -> [
  //   1   5
  //   2   6
  //   3   7
  //   4   8
  // 101 105
  // 102 106
  // 103 107
  // 104 108
  // 201 205
  // 202 206
  // 203 207
  // 204 208];
  //   
  template <typename T>
  inline void transpose_blocks(
    const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & A,
    const size_t k,
    const size_t dim,
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & B);
}

// Implementation
#include <cassert>

template <typename T>
inline void igl::transpose_blocks(
  const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & A,
  const size_t k,
  const size_t dim,
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & B)
{
  // Eigen matrices must be 2d so dim must be only 1 or 2
  assert(dim == 1 || dim == 2);
  // Output is not allowed to be input
  assert(&A != &B);


  // block height, width, and number of blocks
  int m,n;
  if(dim == 1)
  {
    m = A.rows()/k;
    n = A.cols();
  }else// dim == 2
  {
    m = A.rows();
    n = A.cols()/k;
  }

  // resize output
  if(dim == 1)
  {
    B.resize(n*k,m);
  }else//dim ==2
  {
    B.resize(n,m*k);
  }

  // loop over blocks
  for(int b = 0;b<(int)k;b++)
  {
    if(dim == 1)
    {
      B.block(b*n,0,n,m) = A.block(b*m,0,m,n).transpose();
    }else//dim ==2
    {
      B.block(0,b*m,n,m) = A.block(0,b*n,m,n).transpose();
    }
  }
}
#endif
