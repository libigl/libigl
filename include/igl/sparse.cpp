// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "sparse.h"

#include <iostream>
#include <vector>

template <class IndexVector, class ValueVector, typename T>
IGL_INLINE void igl::sparse(
  const IndexVector & I,
  const IndexVector & J,
  const ValueVector & V,
  Eigen::SparseMatrix<T>& X)
{
  size_t m = (size_t)I.maxCoeff()+1;
  size_t n = (size_t)J.maxCoeff()+1;
  return igl::sparse(I,J,V,m,n,X);
}

#include "verbose.h"
template <class IndexVector, class ValueVector, typename T>
IGL_INLINE void igl::sparse(
  const IndexVector & I,
  const IndexVector & J,
  const ValueVector & V,
  const size_t m,
  const size_t n,
  Eigen::SparseMatrix<T>& X)
{
  using namespace std;
  using namespace Eigen;
  assert((int)I.maxCoeff() < (int)m);
  assert((int)I.minCoeff() >= 0);
  assert((int)J.maxCoeff() < (int)n);
  assert((int)J.minCoeff() >= 0);
  assert(I.size() == J.size());
  assert(J.size() == V.size());
  // Really we just need .size() to be the same, but this is safer
  assert(I.rows() == J.rows());
  assert(J.rows() == V.rows());
  assert(I.cols() == J.cols());
  assert(J.cols() == V.cols());
  //// number of values
  //int nv = V.size();

  //Eigen::DynamicSparseMatrix<T, Eigen::RowMajor> dyn_X(m,n);
  //// over estimate the number of entries
  //dyn_X.reserve(I.size());
  //for(int i = 0;i < nv;i++)
  //{
  //  dyn_X.coeffRef((int)I(i),(int)J(i)) += (T)V(i);
  //}
  //X = Eigen::SparseMatrix<T>(dyn_X);
  vector<Triplet<T> > IJV;
  IJV.reserve(I.size());
  for(int x = 0;x<I.size();x++)
  {
    IJV.push_back(Triplet<T >(I(x),J(x),V(x)));
  }
  X.resize(m,n);
  X.setFromTriplets(IJV.begin(),IJV.end());
}

template <typename DerivedD, typename T>
IGL_INLINE void igl::sparse(
  const Eigen::PlainObjectBase<DerivedD>& D,
  Eigen::SparseMatrix<T>& X)
{
  assert(false && "Obsolete. Just call D.sparseView() directly");
  using namespace std;
  using namespace Eigen;
  vector<Triplet<T> > DIJV;
  const int m = D.rows();
  const int n = D.cols();
  for(int i = 0;i<m;i++)
  {
    for(int j = 0;j<n;j++)
    {
      if(D(i,j)!=0)
      {
        DIJV.push_back(Triplet<T>(i,j,D(i,j)));
      }
    }
  }
  X.resize(m,n);
  X.setFromTriplets(DIJV.begin(),DIJV.end());
}

template <typename DerivedD>
IGL_INLINE Eigen::SparseMatrix<typename DerivedD::Scalar > igl::sparse(
  const Eigen::PlainObjectBase<DerivedD>& D)
{
  assert(false && "Obsolete. Just call D.sparseView() directly");
  Eigen::SparseMatrix<typename DerivedD::Scalar > X;
  igl::sparse(D,X);
  return X;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::sparse<Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, double>(Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::Matrix<double, -1, 1, 0, -1, 1> const&, size_t, size_t, Eigen::SparseMatrix<double, 0, int>&);
template void igl::sparse<Eigen::Matrix<double, -1, -1, 0, -1, -1>, double>(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::SparseMatrix<double, 0, int>&);
template Eigen::SparseMatrix<Eigen::Matrix<double, -1, -1, 0, -1, -1>::Scalar, 0, int> igl::sparse<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&);
template Eigen::SparseMatrix<Eigen::Matrix<double, -1, 1, 0, -1, 1>::Scalar, 0, int> igl::sparse<Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&);
template void igl::sparse<Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<double>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >, double>(Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<double>, Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::SparseMatrix<double, 0, int>&);
template void igl::sparse<Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, double>(Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::Matrix<double, -1, 1, 0, -1, 1> const&, Eigen::SparseMatrix<double, 0, int>&);

#ifdef _WIN32
template void igl::sparse<Eigen::Matrix<int,-1,1,0,-1,1>,Eigen::Matrix<double,-1,1,0,-1,1>,std::complex<double> >(Eigen::Matrix<int,-1,1,0,-1,1> const &,Eigen::Matrix<int,-1,1,0,-1,1> const &,Eigen::Matrix<double,-1,1,0,-1,1> const &,unsigned long long int,unsigned long long int,Eigen::SparseMatrix<std::complex<double>,0,int> &);
#endif

#ifndef _WIN32
template void igl::sparse<Eigen::Matrix<int,-1,1,0,-1,1>,Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<double>,Eigen::Matrix<double,-1,1,0,-1,1> >,double>(Eigen::Matrix<int,-1,1,0,-1,1> const&,Eigen::Matrix<int,-1,1,0,-1,1> const&,Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<double>,Eigen::Matrix<double,-1,1,0,-1,1> > const&,unsigned long,unsigned long,Eigen::SparseMatrix<double,0,int>&);
template void igl::sparse<Eigen::Matrix<int,-1,1,0,-1,1>,Eigen::Matrix<double,-1,1,0,-1,1>,std::complex<double> >(Eigen::Matrix<int,-1,1,0,-1,1> const&,Eigen::Matrix<int,-1,1,0,-1,1> const&,Eigen::Matrix<double,-1,1,0,-1,1> const&,unsigned long,unsigned long,Eigen::SparseMatrix<std::complex<double>,0,int>&);
#endif

#endif
