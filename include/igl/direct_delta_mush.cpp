// This file is part of libigl, a simple C++ geometry processing library.
//
// Copyright (C) 2020 Xiangyu Kong <xiangyu.kong@mail.utoronto.ca>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "direct_delta_mush.h"
#include "cotmatrix.h"

#include <iostream>

// TODOs
// 1. U_precomputed, Psi, Omega should be #V by 10 instead of #V by 16
// 2. Vectorize Psi computation.

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedC,
  typename DerivedE,
  typename DerivedOmega,
  typename DerivedU>
IGL_INLINE void igl::direct_delta_mush(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  const Eigen::MatrixBase<DerivedC> & C,
  const Eigen::MatrixBase<DerivedE> & E,
  const std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d> > & T,
  const Eigen::MatrixBase<DerivedOmega> & Omega,
  Eigen::PlainObjectBase<DerivedU> & U)
{
  using namespace std;
  using namespace Eigen;
  cout << "START DDM" << endl;

  cout << "V: " << V.rows() << " x " << V.cols() << " Sum: " << V.sum()
       << "\nF: " << F.rows() << " x " << F.cols() << " Sum: " << F.sum()
       << "\nC: " << C.rows() << " x " << C.cols() << " Sum: " << C.sum()
       << "\nE: " << E.rows() << " x " << E.cols() << " Sum: " << E.sum()
       << "\nT: " << T.size()
       << "\nOmega: " << Omega.rows() << " x " << Omega.cols() << " Sum: " << Omega.sum()
       << endl;

  assert(V.cols() == 3 && "V should contain 3D positions.");
  assert(F.cols() == 3 && "F should contain triangles.");
  assert(C.cols() == 3 && "C should contain 3D bone endpoint positions.");
  assert(E.cols() == 2 && "E should contain 2 endpoint indices forming bone edges.");
  assert(E.rows() == T.size() && "E.rows() should equal to T.size()");
  assert(Omega.rows() == V.rows() && "Omega contain the same number of rows as V.");
  assert(Omega.cols() == T.size() * 16 && "Omega should have #T*10 columns.");

  int n = V.rows();
  int m = T.size();

  Eigen::MatrixXd U_homogeneous(n, 4);
  U_homogeneous << V, Eigen::VectorXd::Ones(n);
  U.resize(n, 3);

  // R matrix
  Eigen::MatrixXd R_matrices(n, 9);
  Eigen::MatrixXd T_matrices(n, 3);
  for (int i = 0; i < n; i++)
  {
    Eigen::MatrixXd Q_mat = Eigen::MatrixXd::Zero(4, 4);
    for (int j = 0; j < m; j++)
    {
      Eigen::MatrixXd Omega_curr = Omega.block(i, j * 16, 1, 16);
      Omega_curr.resize(4, 4);
      Eigen::Affine3d M_curr = T[j];
      Q_mat += M_curr.matrix() * Omega_curr;
    }
    Eigen::MatrixXd Q_i = Q_mat.block(0, 0, 3, 3);
    Eigen::MatrixXd q_i = Q_mat.block(0, 3, 3, 1);
    Eigen::MatrixXd p_i = Q_mat.block(3, 0, 1, 3); // .transpose()

    Eigen::MatrixXd SVD_i = Q_i - q_i * p_i;
    Eigen::JacobiSVD<MatrixXd> svd;
    svd.compute(SVD_i, Eigen::ComputeFullU | Eigen::ComputeFullV );

    // rotation and translation
    Eigen::MatrixXd R_i = svd.matrixU() * svd.matrixV().transpose();
    Eigen::VectorXd t_i = q_i - R_i * p_i.transpose();

    // Gamma
    Eigen::MatrixXd Gamma_i(3, 4);
    Gamma_i.block(0, 0, 3, 3) = R_i;
    Gamma_i.block(0, 3, 3, 1) = t_i;

    // transform
    Eigen::VectorXd u_i = U_homogeneous.row(i);
    U.row(i) = Gamma_i * u_i;
  }

  cout << "END DDM" << endl;
}

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedC,
  typename DerivedE,
  typename DerivedW,
  typename DerivedOmega>
IGL_INLINE void igl::direct_delta_mush_precomputation(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  const Eigen::MatrixBase<DerivedC> & C,
  const Eigen::MatrixBase<DerivedE> & E,
  const Eigen::SparseMatrix<DerivedW> & W,
  const int p,
  const typename DerivedV::Scalar lambda,
  const typename DerivedV::Scalar kappa,
  const typename DerivedV::Scalar alpha,
  Eigen::PlainObjectBase<DerivedOmega> & Omega)
{
  using namespace std;
  using namespace Eigen;

  cout << "START DDM Precomputation" << endl;

  cout << "Using params:"
       << "\np: " << p
       << "\nlambda: " << lambda
       << "\nkappa: " << kappa
       << endl;
  cout << "V: " << V.rows() << " x " << V.cols() << " Sum: " << V.sum()
       << "\nF: " << F.rows() << " x " << F.cols() << " Sum: " << F.sum()
       << "\nC: " << C.rows() << " x " << C.cols() << " Sum: " << C.sum()
       << "\nE: " << E.rows() << " x " << E.cols() << " Sum: " << E.sum()
       << "\nW: " << W.rows() << " x " << W.cols() << " Sum: " << W.sum()
       << endl;

  assert(V.cols() == 3 && "V should contain 3D positions.");
  assert(F.cols() == 3 && "F should contain triangles.");
  assert(C.cols() == 3 && "C should contain 3D bone endpoint positions.");
  assert(E.cols() == 2 && "E should contain 2 endpoint indices forming bone edges.");
  assert(W.rows() == V.rows() && "W.rows() should be equal to V.rows().");
  assert(W.cols() == E.rows() && "W.cols() should be equal to E.rows().");
  assert(p > 0 && "Laplacian iteration p should be positive.");
  assert(lambda > 0 && "lambda should be positive.");
  assert(kappa > 0 && kappa < lambda && "kappa should be positive and less than lambda.");
  assert(alpha >= 0 && alpha < 1 && "alpha should be non-negative and less than 1.");

  const int n = V.rows();
  const int m = E.rows();

  // U: #V by 4, homogeneous version of V
  Eigen::MatrixXd U_homogeneous(n, 4);
  U_homogeneous << V, Eigen::VectorXd::Ones(n);

  // Identity of #V by #V
  Eigen::SparseMatrix<double> I(n, n);
  I.setIdentity();

  // Laplacian: L_bar = L \times D_L^{-1}
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);
  Eigen::MatrixXd D_L = L.diagonal().asDiagonal();
  Eigen::SparseMatrix<double> D_L_sparse = D_L.sparseView();
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
  ldlt.compute(D_L_sparse); // this will take quite a while to compute
  Eigen::SparseMatrix<double> L_bar = ldlt.solve(L).transpose();

  // Implicitly and iteratively solve
  // w'_{ij} = \sum_{k=1}^{n}{C_{ki} w_{kj}}, C = (I + kappa L_bar)^{-p}:
  // W' = C^T \times W
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt_W_prime;
  Eigen::SparseMatrix<double> C_calc((I + kappa * L_bar).transpose());
  Eigen::SparseMatrix<double> W_prime(W);
  ldlt_W_prime.compute(C_calc);
  cout << "computing W" << endl;
  for (int iter = 0; iter < p; iter++)
  {
    W_prime.makeCompressed();
    W_prime = ldlt_W_prime.solve(W_prime);
  }
  cout << "W_prime: " << W_prime.rows() << " x " << W_prime.cols()
       << " Sum: " << W_prime.sum() << endl;

  // Psi was hard to solve iteratively since i couldn't express u_k \times u_k^T as matrix form
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt_U_tilde;
  Eigen::SparseMatrix<double> B_calc((I + lambda * L_bar).transpose());
  ldlt_U_tilde.compute(B_calc);

  Eigen::SparseMatrix<double> U_tilde = U_homogeneous.sparseView();
  cout << "computing U_tilde" << endl;
  for (int i = 0; i < p; i++)
  {
    U_tilde.makeCompressed();
    U_tilde = ldlt_U_tilde.solve(U_tilde);
  }
  cout << "U_tilde: " << U_tilde.rows() << " x " << U_tilde.cols()
       << " Sum: " << U_tilde.sum() << endl;

  // TODO: sparse didn't work here - malloc error
  Eigen::MatrixXd U_tilde_dense = U_tilde.toDense();
  Eigen::MatrixXd B_inv_dense = U_homogeneous.transpose().householderQr().solve(
    U_tilde_dense.transpose());
  Eigen::SparseMatrix<double> B_inv = B_inv_dense.sparseView();
  // Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> qr_B;
  // Eigen::SparseMatrix<double> U_sparse = U_homogeneous.sparseView();
  // Eigen::SparseMatrix<double> V_sparse_transpose(U_sparse.transpose());
  // Eigen::SparseMatrix<double> V_tilde_transpose(U_tilde.transpose());
  // V_sparse_transpose.makeCompressed();
  // qr_B.compute(V_sparse_transpose);
  // Eigen::SparseMatrix<double> B_inv = qr_B.solve(V_tilde_transpose);
  cout << "B_inv: " << B_inv.rows() << " x " << B_inv.cols() << " Sum: " << B_inv.sum() << endl;

  // U_precomputed: #V by 16(10)
  // U_precomputed.row(i) = u_i \dot u_i^T \in R^{4 x 4}
  Eigen::MatrixXd U_precomputed(n, 16);
  for (int k = 0; k < n; k++)
  {
    Eigen::MatrixXd u_full = U_homogeneous.row(k).transpose() * U_homogeneous.row(k);
    Eigen::VectorXd u_full_vector(
      Eigen::Map<Eigen::VectorXd>(u_full.data(), u_full.cols() * u_full.rows()));
    U_precomputed.row(k) = u_full_vector;
  }
  cout << "U_precomputed: " << U_precomputed.rows() << " x " << U_precomputed.cols()
       << " Sum: " << U_precomputed.sum() << endl;

  // Psi: #V by #T*16 (10) of \Psi_{ij}s.
  // this takes a while since it is not vectorized
  // Eigen::MatrixXd Psi = Eigen::MatrixXd::Ones(n, m * 16); // for debugging to skip computation
  Eigen::MatrixXd Psi(n, m * 16);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      Eigen::VectorXd Psi_curr = Eigen::VectorXd::Zero(16);
      for (int k = 0; k < n; k++)
      {
        Psi_curr += B_inv.coeff(k, i) * W.coeff(k, j) * U_precomputed.row(k);
      }
      Psi.block(i, j * 16, 1, 16) = Psi_curr.transpose();
    }
  }
  cout << "Psi: " << Psi.rows() << " x " << Psi.cols() << " Sum: " << Psi.sum() << endl;

  // precomputed P matrix: 4 by 4
  // p_i p_i^T , p_i
  // p_i^T     , 1
  // p_i: sum_{j=1}^{n} Psi_{ij} top right 3 by 1 column
  Eigen::MatrixXd P_vectors(n, 16);
  for (int i = 0; i < n; i++)
  {
    Eigen::Vector3d p_i = Eigen::Vector3d::Zero(3);
    for (int j = 0; j < m; j++)
    {
      Eigen::Vector3d p_i_curr(3);
      p_i_curr << Psi(i, j * 16 + 3), Psi(i, j * 16 + 7), Psi(i, j * 16 + 11);
      p_i += p_i_curr;
    }
    Eigen::MatrixXd p_matrix(4, 4);
    p_matrix.block(0, 0, 3, 3) = p_i * p_i.transpose();
    p_matrix.block(3, 0, 1, 3) = p_i.transpose();
    p_matrix.block(0, 3, 3, 1) = p_i;
    p_matrix(3, 3) = 1;
    P_vectors.row(i) = Eigen::Map<Eigen::VectorXd>(
      p_matrix.data(), p_matrix.cols() * p_matrix.rows());
  }
  cout << "P_vectors: " << P_vectors.rows() << " x " << P_vectors.cols() << " Sum: "
       << P_vectors.sum() << endl;

  // Omega
  Omega.resize(n, m * 16);
  for (int i = 0; i < n; i++)
  {
    Eigen::MatrixXd p_vector = P_vectors.row(i);
    for (int j = 0; j < m; j++)
    {
      Eigen::MatrixXd Omega_curr(1, 16);
      Eigen::MatrixXd Psi_curr = Psi.block(i, j * 16, 1, 16);
      Omega_curr = (1 - alpha) * Psi_curr + alpha * W_prime.coeff(i, j) * p_vector;
      Omega.block(i, j * 16, 1, 16) = Omega_curr;
    }
  }
  cout << "Omega: " << Omega.rows() << " x " << Omega.cols() << " Sum: " << Omega.sum() << endl;

  cout << "END DDM Precomputation" << endl;
}

template <
  typename DerivedOmega,
  typename DerivedU>
IGL_INLINE void igl::direct_delta_mush_pose_evaluation(
  const std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d> > & T,
  const Eigen::MatrixBase<DerivedOmega> & Omega,
  Eigen::PlainObjectBase<DerivedU> & U)
{
  using namespace std;
  using namespace Eigen;

  cout << "START DDM Pose Eval" << endl;
  // not sure what this is
  cout << "END DDM Pose Eval" << endl;
}

#ifdef IGL_STATIC_LIBRARY

// Explicit template instantiation
template void
igl::direct_delta_mush<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(
  Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &,
  Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &,
  Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &,
  Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &,
  std::__1::vector<Eigen::Transform<double, 3, 2, 0>, Eigen::aligned_allocator<Eigen::Transform<double, 3, 2, 0> > > const &,
  Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &,
  Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > &);

template void
igl::direct_delta_mush_precomputation<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(
  Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &,
  Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &,
  Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &,
  Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &,
  Eigen::SparseMatrix<double, 0, int> const &, int,
  Eigen::Matrix<double, -1, -1, 0, -1, -1>::Scalar,
  Eigen::Matrix<double, -1, -1, 0, -1, -1>::Scalar,
  Eigen::Matrix<double, -1, -1, 0, -1, -1>::Scalar,
  Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > &);

template void
igl::direct_delta_mush_pose_evaluation<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(
  std::__1::vector<Eigen::Transform<double, 3, 2, 0>, Eigen::aligned_allocator<Eigen::Transform<double, 3, 2, 0> > > const &,
  Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &,
  Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > &);

#endif