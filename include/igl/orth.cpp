#include "orth.h"

// Broken Implementation
IGL_INLINE void igl::orth(const Eigen::MatrixXd &A, Eigen::MatrixXd &Q)
{

  //perform svd on A = U*S*V' (V is not computed and only the thin U is computed)
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU );
  Eigen::MatrixXd U = svd.matrixU();
  const Eigen::VectorXd S = svd.singularValues();
  
  //get rank of A
  int m = A.rows();
  int n = A.cols();
  double tol = std::max(m,n) * S.maxCoeff() *  2.2204e-16;
  int r = 0;
  for (int i = 0; i < S.rows(); ++r,++i)
  {
    if (S[i] < tol)
      break;
  }
  
  //keep r first columns of U
  Q = U.block(0,0,U.rows(),r);
}
