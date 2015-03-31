// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "slice.h"
#include "colon.h"

#include <vector>

// Bug in unsupported/Eigen/SparseExtra needs iostream first
#include <iostream>
#include <unsupported/Eigen/SparseExtra>

template <typename T>
IGL_INLINE void igl::slice(
  const Eigen::SparseMatrix<T>& X,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & C,
  Eigen::SparseMatrix<T>& Y)
{
  int xm = X.rows();
  int xn = X.cols();
  int ym = R.size();
  int yn = C.size();

  // special case when R or C is empty
  if(ym == 0 || yn == 0)
  {
    Y.resize(ym,yn);
    return;
  }

  assert(R.minCoeff() >= 0);
  assert(R.maxCoeff() < xm);
  assert(C.minCoeff() >= 0);
  assert(C.maxCoeff() < xn);

  // initialize row and col permutation vectors
  Eigen::VectorXi rowIndexVec = Eigen::VectorXi::LinSpaced(xm,0,xm);
  Eigen::VectorXi rowPermVec  = Eigen::VectorXi::LinSpaced(xm,0,xm);
  for(int i=0;i<ym;i++)
  {
    int pos = rowIndexVec.coeffRef(R(i));
    if(pos != i)
    {
      int& val = rowPermVec.coeffRef(i);
      std::swap(rowIndexVec.coeffRef(val),rowIndexVec.coeffRef(R(i)));
      std::swap(rowPermVec.coeffRef(i),rowPermVec.coeffRef(pos));
    }
  }
  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic,int> rowPerm(rowIndexVec);

  Eigen::VectorXi colIndexVec = Eigen::VectorXi::LinSpaced(xn,0,xn);
  Eigen::VectorXi colPermVec = Eigen::VectorXi::LinSpaced(xn,0,xn);
  for(int i=0;i<yn;i++)
  {
    int pos = colIndexVec.coeffRef(C(i));
    if(pos != i)
    {
      int& val = colPermVec.coeffRef(i);
      std::swap(colIndexVec.coeffRef(val),colIndexVec.coeffRef(C(i)));
      std::swap(colPermVec.coeffRef(i),colPermVec.coeffRef(pos));
    }
  }
  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic,int> colPerm(colPermVec);

  Eigen::SparseMatrix<T> M = (rowPerm * X);
  Y = (M * colPerm).block(0,0,ym,yn);
}

template <typename Mat>
IGL_INLINE void igl::slice(
  const Mat& X,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
  const int dim,
  Mat& Y)
{
  Eigen::VectorXi C;
  switch(dim)
  {
    case 1:
      // boring base case
      if(X.cols() == 0)
      {
        Y.resize(R.size(),0);
        return;
      }
      igl::colon(0,X.cols()-1,C);
      return slice(X,R,C,Y);
    case 2:
      // boring base case
      if(X.rows() == 0)
      {
        Y.resize(0,R.size());
        return;
      }
      igl::colon(0,X.rows()-1,C);
      return slice(X,C,R,Y);
    default:
      assert(false && "Unsupported dimension");
      return;
  }
}

template <typename DerivedX>
IGL_INLINE void igl::slice(
  const Eigen::PlainObjectBase<DerivedX> & X,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & C,
  Eigen::PlainObjectBase<DerivedX> & Y)
{
#ifndef NDEBUG
  int xm = X.rows();
  int xn = X.cols();
#endif
  int ym = R.size();
  int yn = C.size();

  // special case when R or C is empty
  if(ym == 0 || yn == 0)
  {
    Y.resize(ym,yn);
    return;
  }

  assert(R.minCoeff() >= 0);
  assert(R.maxCoeff() < xm);
  assert(C.minCoeff() >= 0);
  assert(C.maxCoeff() < xn);

  // Resize output
  Y.resize(ym,yn);
  // loop over output rows, then columns
  for(int i = 0;i<ym;i++)
  {
    for(int j = 0;j<yn;j++)
    {
      Y(i,j) = X(R(i),C(j));
    }
  }
}

template <typename DerivedX>
IGL_INLINE void igl::slice_mask(
  const Eigen::PlainObjectBase<DerivedX> & X,
  const Eigen::Array<bool,Eigen::Dynamic,1> & R,
  const Eigen::Array<bool,Eigen::Dynamic,1> & C,
  Eigen::PlainObjectBase<DerivedX> & Y)
{
  int xm = X.rows();
  int xn = X.cols();
  int ym = R.count();
  int yn = C.count();
  assert(R.size() == X.rows() && "R.size() should match X.rows()");
  assert(C.size() == X.cols() && "C.size() should match X.cols()");
  Y.resize(ym,yn);
  {
    int yi = 0;
    for(int i = 0;i<xm;i++)
    {
      if(R(i))
      {
        int yj = 0;
        for(int j = 0;j<xn;j++)
        {
          if(C(j))
          {
            Y(yi,yj) = X(i,j);
            yj++;
          }
        }
        yi++;
      }
    }
  }
}

template <typename DerivedX>
IGL_INLINE void igl::slice_mask(
  const Eigen::PlainObjectBase<DerivedX> & X,
  const Eigen::Array<bool,Eigen::Dynamic,1> & R,
  const int dim,
  Eigen::PlainObjectBase<DerivedX> & Y)
{
  int xm = X.rows();
  int xn = X.cols();
  switch(dim)
  {
    case 1:
    {
      const int ym = R.count();
      Y.resize(ym,X.cols());
      assert(X.rows() == R.size() && "X.rows() should match R.size()");
      {
        int yi = 0;
        for(int i = 0;i<X.rows();i++)
        {
          if(R(i))
          {
            Y.row(yi++) = X.row(i);
          }
        }
      }
      return;
    }
    case 2:
    {
      const auto & C = R;
      const int yn = C.count();
      Y.resize(X.rows(),yn);
      assert(X.cols() == R.size() && "X.cols() should match R.size()");
      {
        int yj = 0;
        for(int j = 0;j<X.cols();j++)
        {
          if(C(j))
          {
            Y.col(yj++) = X.col(j);
          }
        }
      }
      return;
    }
    default:
      assert(false && "Unsupported dimension");
      return;
  }
}

template <typename DerivedX>
IGL_INLINE void igl::slice(
  const Eigen::PlainObjectBase<DerivedX> & X,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
  Eigen::PlainObjectBase<DerivedX> & Y)
{
  // phony column indices
  Eigen::Matrix<int,Eigen::Dynamic,1> C;
  C.resize(1);
  C(0) = 0;
  return igl::slice(X,R,C,Y);
}

template <typename DerivedX>
IGL_INLINE Eigen::PlainObjectBase<DerivedX> igl::slice(
  const Eigen::PlainObjectBase<DerivedX> & X,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & R)
{
  Eigen::PlainObjectBase<DerivedX> Y;
  igl::slice(X,R,Y);
  return Y;
}

template <typename DerivedX>
IGL_INLINE Eigen::PlainObjectBase<DerivedX> igl::slice(
  const Eigen::PlainObjectBase<DerivedX>& X,
  const Eigen::Matrix<int,Eigen::Dynamic,1> & R,
  const int dim)
{
  Eigen::PlainObjectBase<DerivedX> Y;
  igl::slice(X,R,dim,Y);
  return Y;
}

template <typename DerivedX>
IGL_INLINE Eigen::PlainObjectBase<DerivedX> igl::slice_mask(
  const Eigen::PlainObjectBase<DerivedX> & X,
  const Eigen::Array<bool,Eigen::Dynamic,1> & R,
  const Eigen::Array<bool,Eigen::Dynamic,1> & C)
{
  Eigen::PlainObjectBase<DerivedX> Y;
  igl::slice_mask(X,R,C,Y);
  return Y;
}

template <typename DerivedX>
IGL_INLINE Eigen::PlainObjectBase<DerivedX> igl::slice_mask(
  const Eigen::PlainObjectBase<DerivedX>& X,
  const Eigen::Array<bool,Eigen::Dynamic,1> & R,
  const int dim)
{
  Eigen::PlainObjectBase<DerivedX> Y;
  igl::slice_mask(X,R,dim,Y);
  return Y;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
// generated by autoexplicit.sh
template void igl::slice<Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
// generated by autoexplicit.sh
template void igl::slice<Eigen::Matrix<float, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<float, -1, 1, 0, -1, 1> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, 1, 0, -1, 1> >&);
// generated by autoexplicit.sh
template void igl::slice<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
// generated by autoexplicit.sh
template void igl::slice<Eigen::Matrix<float, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> >&);
// generated by autoexplicit.sh
template void igl::slice<double>(Eigen::SparseMatrix<double, 0, int> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::SparseMatrix<double, 0, int>&);
template void igl::slice<Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
template void igl::slice<Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::slice<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, int, Eigen::Matrix<double, -1, -1, 0, -1, -1>&);
template void igl::slice<Eigen::SparseMatrix<double, 0, int> >(Eigen::SparseMatrix<double, 0, int> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, int, Eigen::SparseMatrix<double, 0, int>&);
template Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > igl::slice<Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, int);
template Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > igl::slice<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, int);
template void igl::slice<std::complex<double> >(Eigen::SparseMatrix<std::complex<double>, 0, int> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::SparseMatrix<std::complex<double>, 0, int>&);
template void igl::slice<Eigen::Matrix<double, -1, 3, 0, -1, 3> >(Eigen::Matrix<double, -1, 3, 0, -1, 3> const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, int, Eigen::Matrix<double, -1, 3, 0, -1, 3>&);
template Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > igl::slice<Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&);
template Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > igl::slice<Eigen::Matrix<double, -1, 3, 0, -1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, int);
template Eigen::PlainObjectBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > igl::slice<Eigen::Matrix<double, 1, -1, 1, 1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, int);
template Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > igl::slice<Eigen::Matrix<int, -1, 3, 0, -1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, int);
#endif
