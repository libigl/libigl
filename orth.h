//
//  orth.h
//
//  Created by Olga Diamanti on 9/11/11.
//  Copyright (c) 2011 ETH Zurich. All rights reserved.
//

#ifndef IGL_ORTH_H
#define IGL_ORTH_H

#include <Eigen/Core>

namespace igl
{
//  ORTH   Orthogonalization.
//     ORTH(A,Q) produces Q as an orthonormal basis for the range of A.
//     That is, Q'*Q = I, the columns of Q span the same space as 
//     the columns of A, and the number of columns of Q is the 
//     rank of A.
//  
//  
//   The algorithm  uses singular value decomposition, SVD, instead of orthogonal
//   factorization, QR.  This doubles the computation time, but
//   provides more reliable and consistent rank determination.
//   Closely follows MATLAB implementation in orth.m
//  
  inline void orth(const Eigen::MatrixXd &A, Eigen::MatrixXd &Q);
}

// Implementation
inline void igl::orth(const Eigen::MatrixXd &A, Eigen::MatrixXd &Q)
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

#endif
