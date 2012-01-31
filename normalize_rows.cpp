#include "normalize_rows.h"

IGL_INLINE void igl::normalize_rows(const Eigen::MatrixXd & A, Eigen::MatrixXd & B)
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
