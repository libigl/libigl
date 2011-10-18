#ifndef IGL_NORMALIZE_ROWS_H
#define IGL_NORMALIZE_ROWS_H
#include <Eigen/Core>

namespace igl
{
  // Normalize the rows in A so that their lengths are each 1 and place the new
  // entries in B
  // Inputs:
  //   A  #rows by k input matrix
  // Outputs:
  //   B  #rows by k input matrix, can be the same as A
  inline void normalize_rows(const Eigen::MatrixXd & A, Eigen::MatrixXd & B);
}

inline void igl::normalize_rows(const Eigen::MatrixXd & A, Eigen::MatrixXd & B)
{
  // Resize output
  B.resize(A.rows(),A.cols());

  // loop over rows
  for(int i = 0; i < A.rows();i++)
  {
    double length = 0;
    // loop over cols
    for(int j = 0; j < A.cols();j++)
    {
      length += A(i,j)*A(i,j);
    }
    length = sqrt(length);
    // loop over cols
    for(int j = 0; j < A.cols();j++)
    {
      B(i,j) = A(i,j) / length;
    }
  }
}

#endif
