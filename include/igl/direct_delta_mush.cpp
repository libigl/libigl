// This file is part of libigl, a simple C++ geometry processing library.
//
// Copyright (C) 2020 Xiangyu Kong <xiangyu.kong@mail.utoronto.ca>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "direct_delta_mush.h"
#include "cotmatrix.h"
#include "diag.h"


// ===== DEBUG: START
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <iostream>
#include "Timer.h"

using namespace std;
// ===== DEBUG: END

// TODOs:

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedC,
  typename DerivedE,
  typename DerivedW,
  typename DerivedT,
  typename DerivedTAlloc,
  typename DerivedU>
IGL_INLINE void igl::direct_delta_mush(
  const Eigen::MatrixBase<DerivedV> &V,
  const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedC> &C,
  const Eigen::MatrixBase<DerivedE> &E,
  const Eigen::SparseMatrix<DerivedW> &W,
  const std::vector<DerivedT, DerivedTAlloc> &T,
  Eigen::PlainObjectBase<DerivedU> &U)
{
  cout << "START DDM" << endl;
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
  const Eigen::MatrixBase<DerivedV> &V,
  const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedC> &C,
  const Eigen::MatrixBase<DerivedE> &E,
  const Eigen::SparseMatrix<DerivedW> &W,
  const int p,
  const typename DerivedV::Scalar lambda,
  const typename DerivedV::Scalar kappa,
  Eigen::PlainObjectBase<DerivedOmega> &Omega)
{
  assert(kappa < lambda &&
    "kappa needs to be smaller than lambda so that optimization for R_i is well defined");
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

  const int n = V.rows();
  const int m = C.rows();

  // Identity of #V by #V
  Eigen::SparseMatrix<double> I(n, n);
  I.setIdentity();
  cout << "I: " << I.rows() << " x " << I.cols() << " Sum: " << I.sum() << endl;

  // Laplacian: L_bar = L \times D_L^{-1}
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);
  cout << "L: " << L.rows() << " x " << L.cols() << " Sum: " << L.sum() << endl;
  Eigen::MatrixXd D_L = L.diagonal().asDiagonal();
  cout << "D_L: " << D_L.rows() << " x " << D_L.cols() << " Sum: " << D_L.sum() << endl;
  Eigen::SparseMatrix<double> D_L_sparse = D_L.sparseView();
  cout << "D_L_sparse: " << D_L_sparse.rows() << " x " << D_L_sparse.cols() << " Sum: " << D_L_sparse.sum() << endl;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
  ldlt.compute(D_L_sparse);
  Eigen::SparseMatrix<double> L_bar = ldlt.solve(L).transpose();
  cout << "L_bar: " << L_bar.rows() << " x " << L_bar.cols() << " Sum: " << L_bar.sum() << endl;

  // Implicitly and iteratively solve
  // w'_{ij} = \sum_{k=1}^{n}{C_{ki} w_{kj}}, C = (I + kappa L_bar)^{-p}:
  // W' = C^T \times W
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt_C;
  Eigen::SparseMatrix<double> C_calc((I + kappa * L_bar).transpose());
  cout << "C_calc: " << C_calc.rows() << " x " << C_calc.cols() << " Sum: " << C_calc.sum() << endl;
  Eigen::SparseMatrix<double> W_prime(W), W_next;
  cout << "W_prime: " << W_prime.rows() << " x " << W_prime.cols() << " Sum: " << W_prime.sum() << endl;
  ldlt_C.compute(C_calc);
  cout << "computing W'" << endl;
  for (int iter = 0; iter < p; iter++) {
    cout << "iter:" << iter << endl;
    W_prime.makeCompressed();
    W_prime = ldlt_C.solve(W_prime);
  }
  cout << "W_prime: " << W_prime.rows() << " x " << W_prime.cols() << " Sum: " << W_prime.sum() << endl;

  cout << "END DDM Precomputation" << endl;
}

template <
  typename DerivedT,
  typename DerivedTAlloc,
  typename DerivedOmega,
  typename DerivedU>
IGL_INLINE void igl::direct_delta_mush_pose_evaluation(
  const std::vector<DerivedT, DerivedTAlloc> &T,
  const Eigen::MatrixBase<DerivedOmega> &Omega,
  Eigen::PlainObjectBase<DerivedU> &U)
{
  cout << "START DDM Pose Eval" << endl;
  cout << "END DDM Pose Eval" << endl;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::direct_delta_mush<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double, Eigen::Transform<double, 3, 2, 0>, Eigen::aligned_allocator<Eigen::Transform<double, 3, 2, 0> >, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::SparseMatrix<double, 0, int> const&, std::__1::vector<Eigen::Transform<double, 3, 2, 0>, Eigen::aligned_allocator<Eigen::Transform<double, 3, 2, 0> > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::direct_delta_mush_precomputation<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::SparseMatrix<double, 0, int> const&, int, Eigen::Matrix<double, -1, -1, 0, -1, -1>::Scalar, Eigen::Matrix<double, -1, -1, 0, -1, -1>::Scalar, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::direct_delta_mush_pose_evaluation<Eigen::Transform<double, 3, 2, 0>, Eigen::aligned_allocator<Eigen::Transform<double, 3, 2, 0> >, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::__1::vector<Eigen::Transform<double, 3, 2, 0>, Eigen::aligned_allocator<Eigen::Transform<double, 3, 2, 0> > > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif