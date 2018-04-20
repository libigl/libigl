// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "harmonic.h"
#include "adjacency_matrix.h"
#include "cotmatrix.h"
#include "diag.h"
#include "invert_diag.h"
#include "isdiag.h"
#include "massmatrix.h"
#include "min_quad_with_fixed.h"
#include "speye.h"
#include "sum.h"
#include <Eigen/Sparse>

template <
  typename DerivedV,
  typename DerivedF,
  typename Derivedb,
  typename Derivedbc,
  typename DerivedW>
IGL_INLINE bool igl::harmonic(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  const Eigen::MatrixBase<Derivedb> & b,
  const Eigen::MatrixBase<Derivedbc> & bc,
  const int k,
  Eigen::PlainObjectBase<DerivedW> & W)
{
  using namespace Eigen;
  typedef typename DerivedV::Scalar Scalar;
  SparseMatrix<Scalar> L,M;
  cotmatrix(V,F,L);
  if(k>1)
  {
    massmatrix(V,F,MASSMATRIX_TYPE_DEFAULT,M);
  }
  return harmonic(L,M,b,bc,k,W);
}

template <
  typename DerivedF,
  typename Derivedb,
  typename Derivedbc,
  typename DerivedW>
IGL_INLINE bool igl::harmonic(
  const Eigen::MatrixBase<DerivedF> & F,
  const Eigen::MatrixBase<Derivedb> & b,
  const Eigen::MatrixBase<Derivedbc> & bc,
  const int k,
  Eigen::PlainObjectBase<DerivedW> & W)
{
  using namespace Eigen;
  typedef typename Derivedbc::Scalar Scalar;
  SparseMatrix<Scalar> A;
  adjacency_matrix(F,A);
  // sum each row
  SparseVector<Scalar> Asum;
  sum(A,1,Asum);
  // Convert row sums into diagonal of sparse matrix
  SparseMatrix<Scalar> Adiag;
  diag(Asum,Adiag);
  SparseMatrix<Scalar> L = A-Adiag;
  SparseMatrix<Scalar> M;
  speye(L.rows(),M);
  return harmonic(L,M,b,bc,k,W);
}

template <
  typename DerivedF,
  typename Derivedb,
  typename Derivedbc,
  typename VectorIndex,
  typename DerivedW>
IGL_INLINE bool igl::harmonic(
  const Eigen::MatrixBase<DerivedF> & F,
  const Eigen::MatrixBase<Derivedb> & b,
  const Eigen::MatrixBase<Derivedbc> & bc,
  const std::vector<VectorIndex> & holes,
  const int k,
  Eigen::PlainObjectBase<DerivedW> & W)
{
  int n_filled_faces = 0;
  int num_holes = holes.size();
  int real_F_num = F.rows();
  const int V_rows = F.maxCoeff()+1;

  for (int i = 0; i < num_holes; i++)
    n_filled_faces += holes[i].size();
  DerivedF F_filled(n_filled_faces + real_F_num, 3);
  F_filled.topRows(real_F_num) = F;

  int new_vert_id = V_rows;
  int new_face_id = real_F_num;

  for (int i = 0; i < num_holes; i++, new_vert_id++)
  {
    int cur_bnd_size = holes[i].size();
    int it = 0;
    int back = holes[i].size() - 1;
    F_filled.row(new_face_id++) << holes[i][it], holes[i][back], new_vert_id;
    while (it != back)
    {
      F_filled.row(new_face_id++)
          << holes[i][(it + 1)],
          holes[i][(it)], new_vert_id;
      it++;
    }
  }
  assert(new_face_id == F_filled.rows());
  assert(new_vert_id == V_rows + num_holes);

  bool flag = harmonic(F_filled, b, bc, k, W);
  W.conservativeResize(V_rows, 2);
  return flag;
}

template <
  typename DerivedL,
  typename DerivedM,
  typename Derivedb,
  typename Derivedbc,
  typename DerivedW>
IGL_INLINE bool igl::harmonic(
  const Eigen::SparseMatrix<DerivedL> & L,
  const Eigen::SparseMatrix<DerivedM> & M,
  const Eigen::MatrixBase<Derivedb> & b,
  const Eigen::MatrixBase<Derivedbc> & bc,
  const int k,
  Eigen::PlainObjectBase<DerivedW> & W)
{
  const int n = L.rows();
  assert(n == L.cols() && "L must be square");
  assert((k==1 || n == M.cols() ) && "M must be same size as L");
  assert((k==1 || n == M.rows() ) && "M must be square");
  assert((k==1 || igl::isdiag(M))  && "Mass matrix should be diagonal");

  Eigen::SparseMatrix<DerivedL> Q;
  igl::harmonic(L,M,k,Q);

  typedef DerivedL Scalar;
  min_quad_with_fixed_data<Scalar> data;
  min_quad_with_fixed_precompute(Q,b,Eigen::SparseMatrix<Scalar>(),true,data);
  W.resize(n,bc.cols());
  typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> VectorXS;
  const VectorXS B = VectorXS::Zero(n,1);
  for(int w = 0;w<bc.cols();w++)
  {
    const VectorXS bcw = bc.col(w);
    VectorXS Ww;
    if(!min_quad_with_fixed_solve(data,B,bcw,VectorXS(),Ww))
    {
      return false;
    }
    W.col(w) = Ww;
  }
  return true;
}

template <
  typename DerivedL,
  typename DerivedM,
  typename DerivedQ>
IGL_INLINE void igl::harmonic(
  const Eigen::SparseMatrix<DerivedL> & L,
  const Eigen::SparseMatrix<DerivedM> & M,
  const int k,
  Eigen::SparseMatrix<DerivedQ> & Q)
{
  assert(L.rows() == L.cols()&&"L should be square");
  Q = -L;
  if(k == 1) return;
  assert(L.rows() == M.rows()&&"L should match M's dimensions");
  assert(M.rows() == M.cols()&&"M should be square");
  Eigen::SparseMatrix<DerivedM> Mi;
  invert_diag(M,Mi);
  // This is **not** robust for k>2. See KKT system in [Jacobson et al. 2010]
  // of the kharmonic function in gptoolbox
  for(int p = 1;p<k;p++)
  {
    Q = (Q*Mi*-L).eval();
  }
}

#include "find.h"
template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedQ>
IGL_INLINE void igl::harmonic(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  const int k,
  Eigen::SparseMatrix<DerivedQ> & Q)
{
  Eigen::SparseMatrix<DerivedQ> L,M;
  cotmatrix(V,F,L);
  if(k>1)
  {
    massmatrix(V,F,MASSMATRIX_TYPE_DEFAULT,M);
  }
  return harmonic(L,M,k,Q);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template bool igl::harmonic<Eigen::Matrix<double, -1, -1, 1, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
// generated by autoexplicit.sh
template bool igl::harmonic<Eigen::Matrix<double, -1, -1, 1, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 1, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 1, -1, -1> >&);
// generated by autoexplicit.sh
template bool igl::harmonic<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
// generated by autoexplicit.sh
template void igl::harmonic<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double>(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, int, Eigen::SparseMatrix<double, 0, int>&);
// generated by autoexplicit.sh
template bool igl::harmonic<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
// generated by autoexplicit.sh
template bool igl::harmonic<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template bool igl::harmonic<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
