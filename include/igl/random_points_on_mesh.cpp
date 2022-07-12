// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "random_points_on_mesh.h"
#include "doublearea.h"
#include "cumsum.h"
#include "histc.h"
#include <iostream>
#include <cassert>
#include <random>

template <typename DerivedV, typename DerivedF, typename DerivedB, typename DerivedFI>
IGL_INLINE void igl::random_points_on_mesh(
  const int n,
  const Eigen::MatrixBase<DerivedV > & V,
  const Eigen::MatrixBase<DerivedF > & F,
  Eigen::PlainObjectBase<DerivedB > & B,
  Eigen::PlainObjectBase<DerivedFI > & FI)
{
  return igl::random_points_on_mesh(n,V,F,std::minstd_rand(std::rand()),B,FI);
}

template <typename URBG, typename DerivedV, typename DerivedF, typename DerivedB, typename DerivedFI>
IGL_INLINE void igl::random_points_on_mesh(
  URBG && urbg,
  const int n,
  const Eigen::MatrixBase<DerivedV > & V,
  const Eigen::MatrixBase<DerivedF > & F,
  Eigen::PlainObjectBase<DerivedB > & B,
  Eigen::PlainObjectBase<DerivedFI > & FI)
{
  using namespace Eigen;
  using namespace std;
  typedef typename DerivedV::Scalar Scalar;
  typedef Matrix<Scalar,Dynamic,1> VectorXs;
  VectorXs A;
  doublearea(V,F,A);
  // Should be traingle mesh. Although Turk's method 1 generalizes...
  assert(F.cols() == 3);
  VectorXs C;
  VectorXs A0(A.size()+1);
  A0(0) = 0;
  A0.bottomRightCorner(A.size(),1) = A;
  // Even faster would be to use the "Alias Table Method"
  cumsum(A0,1,C);
  const Scalar Cmax = C(C.size()-1);
  assert(Cmax > 0 && "Total surface area should be positive");
  // Why is this more accurate than `C /= C(C.size()-1)` ?
  for(int i = 0;i<C.size();i++) { C(i) = C(i)/Cmax; }
  std::uniform_real_distribution<float> dis(-1.0, 1.0);
  const VectorXs R = (VectorXs::NullaryExpr(n,1,[&](){return dis(urbg);}).array() + 1.)/2.;
  assert(R.minCoeff() >= 0);
  assert(R.maxCoeff() <= 1);
  histc(R,C,FI);
  FI = FI.array().min(F.rows() - 1); // fix the bin when R(i) == 1 exactly
  const VectorXs S = (VectorXs::NullaryExpr(n,1,[&](){return dis(urbg);}).array() + 1.)/2.;
  const VectorXs T = (VectorXs::NullaryExpr(n,1,[&](){return dis(urbg);}).array() + 1.)/2.;
  B.resize(n,3);
  B.col(0) = 1.-T.array().sqrt();
  B.col(1) = (1.-S.array()) * T.array().sqrt();
  B.col(2) = S.array() * T.array().sqrt();
}

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedB,
  typename DerivedFI,
  typename DerivedX>
IGL_INLINE void igl::random_points_on_mesh(
  const int n,
  const Eigen::MatrixBase<DerivedV > & V,
  const Eigen::MatrixBase<DerivedF > & F,
  Eigen::PlainObjectBase<DerivedB > & B,
  Eigen::PlainObjectBase<DerivedFI > & FI,
  Eigen::PlainObjectBase<DerivedX> & X)
{
  return igl::random_points_on_mesh(std::minstd_rand(std::rand()),n,V,F,B,FI,X);
}

template <
  typename URBG,
  typename DerivedV,
  typename DerivedF,
  typename DerivedB,
  typename DerivedFI,
  typename DerivedX>
IGL_INLINE void igl::random_points_on_mesh(
  URBG && urbg,
  const int n,
  const Eigen::MatrixBase<DerivedV > & V,
  const Eigen::MatrixBase<DerivedF > & F,
  Eigen::PlainObjectBase<DerivedB > & B,
  Eigen::PlainObjectBase<DerivedFI > & FI,
  Eigen::PlainObjectBase<DerivedX> & X)
{
  random_points_on_mesh(urbg,n,V,F,B,FI);
  X = DerivedX::Zero(B.rows(),V.cols());
  for(int x = 0;x<B.rows();x++)
  {
    for(int b = 0;b<B.cols();b++)
    {
        auto fi = FI(x);
        auto f = F(fi, b);
        auto bval = B(x, b);
        X.row(x) += bval*V.row(f);
    }
  }
}

template <typename DerivedV, typename DerivedF, typename ScalarB, typename DerivedFI>
IGL_INLINE void igl::random_points_on_mesh(
  const int n,
  const Eigen::MatrixBase<DerivedV > & V,
  const Eigen::MatrixBase<DerivedF > & F,
  Eigen::SparseMatrix<ScalarB > & B,
  Eigen::PlainObjectBase<DerivedFI > & FI)
{
  return igl::random_points_on_mesh(std::minstd_rand(std::rand()),n,V,F,B,FI);
}

template <typename URBG, typename DerivedV, typename DerivedF, typename ScalarB, typename DerivedFI>
IGL_INLINE void igl::random_points_on_mesh(
  URBG && urbg,
  const int n,
  const Eigen::MatrixBase<DerivedV > & V,
  const Eigen::MatrixBase<DerivedF > & F,
  Eigen::SparseMatrix<ScalarB > & B,
  Eigen::PlainObjectBase<DerivedFI > & FI)
{
  using namespace Eigen;
  using namespace std;
  Matrix<ScalarB,Dynamic,3> BC;
  random_points_on_mesh(urbg,n,V,F,BC,FI);
  vector<Triplet<ScalarB> > BIJV;
  BIJV.reserve(n*3);
  for(int s = 0;s<n;s++)
  {
    for(int c = 0;c<3;c++)
    {
      assert(FI(s) < F.rows());
      assert(FI(s) >= 0);
      const int v = F(FI(s),c);
      BIJV.push_back(Triplet<ScalarB>(s,v,BC(s,c)));
    }
  }
  B.resize(n,V.rows());
  B.reserve(n*3);
  B.setFromTriplets(BIJV.begin(),BIJV.end());
}

#ifdef IGL_STATIC_LIBRARY
template void igl::random_points_on_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::random_points_on_mesh<std::shuffle_order_engine<std::linear_congruential_engine<unsigned long, 16807ul, 0ul, 2147483647ul>, 256ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::shuffle_order_engine<std::linear_congruential_engine<unsigned long, 16807ul, 0ul, 2147483647ul>, 256ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::random_points_on_mesh<std::discard_block_engine<std::subtract_with_carry_engine<unsigned long, 48ul, 5ul, 12ul>, 389ul, 11ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::discard_block_engine<std::subtract_with_carry_engine<unsigned long, 48ul, 5ul, 12ul>, 389ul, 11ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::random_points_on_mesh<std::discard_block_engine<std::subtract_with_carry_engine<unsigned long, 24ul, 10ul, 24ul>, 223ul, 23ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::discard_block_engine<std::subtract_with_carry_engine<unsigned long, 24ul, 10ul, 24ul>, 223ul, 23ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::random_points_on_mesh<std::subtract_with_carry_engine<unsigned long, 48ul, 5ul, 12ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::subtract_with_carry_engine<unsigned long, 48ul, 5ul, 12ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::random_points_on_mesh<std::subtract_with_carry_engine<unsigned long, 24ul, 10ul, 24ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::subtract_with_carry_engine<unsigned long, 24ul, 10ul, 24ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::random_points_on_mesh<std::mersenne_twister_engine<unsigned long, 64ul, 312ul, 156ul, 31ul, 13043109905998158313ul, 29ul, 6148914691236517205ul, 17ul, 8202884508482404352ul, 37ul, 18444473444759240704ul, 43ul, 6364136223846793005ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::mersenne_twister_engine<unsigned long, 64ul, 312ul, 156ul, 31ul, 13043109905998158313ul, 29ul, 6148914691236517205ul, 17ul, 8202884508482404352ul, 37ul, 18444473444759240704ul, 43ul, 6364136223846793005ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::random_points_on_mesh<std::mersenne_twister_engine<unsigned long, 32ul, 624ul, 397ul, 31ul, 2567483615ul, 11ul, 4294967295ul, 7ul, 2636928640ul, 15ul, 4022730752ul, 18ul, 1812433253ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::mersenne_twister_engine<unsigned long, 32ul, 624ul, 397ul, 31ul, 2567483615ul, 11ul, 4294967295ul, 7ul, 2636928640ul, 15ul, 4022730752ul, 18ul, 1812433253ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::random_points_on_mesh<std::linear_congruential_engine<unsigned long, 48271ul, 0ul, 2147483647ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::linear_congruential_engine<unsigned long, 48271ul, 0ul, 2147483647ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::random_points_on_mesh<std::linear_congruential_engine<unsigned long, 16807ul, 0ul, 2147483647ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::linear_congruential_engine<unsigned long, 16807ul, 0ul, 2147483647ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::random_points_on_mesh<std::linear_congruential_engine<unsigned long, 48271ul, 0ul, 2147483647ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 3, 1, -1, 3> >(std::linear_congruential_engine<unsigned long, 48271ul, 0ul, 2147483647ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&);
template void igl::random_points_on_mesh<std::shuffle_order_engine<std::linear_congruential_engine<unsigned long, 16807ul, 0ul, 2147483647ul>, 256ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 3, 1, -1, 3> >(std::shuffle_order_engine<std::linear_congruential_engine<unsigned long, 16807ul, 0ul, 2147483647ul>, 256ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&);
template void igl::random_points_on_mesh<std::linear_congruential_engine<unsigned long, 16807ul, 0ul, 2147483647ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 3, 1, -1, 3> >(std::linear_congruential_engine<unsigned long, 16807ul, 0ul, 2147483647ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&);
template void igl::random_points_on_mesh<std::subtract_with_carry_engine<unsigned long, 48ul, 5ul, 12ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 3, 1, -1, 3> >(std::subtract_with_carry_engine<unsigned long, 48ul, 5ul, 12ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&);
template void igl::random_points_on_mesh<std::discard_block_engine<std::subtract_with_carry_engine<unsigned long, 48ul, 5ul, 12ul>, 389ul, 11ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 3, 1, -1, 3> >(std::discard_block_engine<std::subtract_with_carry_engine<unsigned long, 48ul, 5ul, 12ul>, 389ul, 11ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&);
template void igl::random_points_on_mesh<std::subtract_with_carry_engine<unsigned long, 24ul, 10ul, 24ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 3, 1, -1, 3> >(std::subtract_with_carry_engine<unsigned long, 24ul, 10ul, 24ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&);
template void igl::random_points_on_mesh<std::discard_block_engine<std::subtract_with_carry_engine<unsigned long, 24ul, 10ul, 24ul>, 223ul, 23ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 3, 1, -1, 3> >(std::discard_block_engine<std::subtract_with_carry_engine<unsigned long, 24ul, 10ul, 24ul>, 223ul, 23ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&);
template void igl::random_points_on_mesh<std::mersenne_twister_engine<unsigned long, 32ul, 624ul, 397ul, 31ul, 2567483615ul, 11ul, 4294967295ul, 7ul, 2636928640ul, 15ul, 4022730752ul, 18ul, 1812433253ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 3, 1, -1, 3> >(std::mersenne_twister_engine<unsigned long, 32ul, 624ul, 397ul, 31ul, 2567483615ul, 11ul, 4294967295ul, 7ul, 2636928640ul, 15ul, 4022730752ul, 18ul, 1812433253ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&);
template void igl::random_points_on_mesh<std::mersenne_twister_engine<unsigned long, 64ul, 312ul, 156ul, 31ul, 13043109905998158313ul, 29ul, 6148914691236517205ul, 17ul, 8202884508482404352ul, 37ul, 18444473444759240704ul, 43ul, 6364136223846793005ul>&, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 3, 1, -1, 3> >(std::mersenne_twister_engine<unsigned long, 64ul, 312ul, 156ul, 31ul, 13043109905998158313ul, 29ul, 6148914691236517205ul, 17ul, 8202884508482404352ul, 37ul, 18444473444759240704ul, 43ul, 6364136223846793005ul>&, int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&);
template void igl::random_points_on_mesh<std::linear_congruential_engine<unsigned long, 48271ul, 0ul, 2147483647ul>&, Eigen::Matrix<float, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<float, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<float, -1, 3, 1, -1, 3> >(std::linear_congruential_engine<unsigned long, 48271ul, 0ul, 2147483647ul>&, int, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> >&);
#endif
