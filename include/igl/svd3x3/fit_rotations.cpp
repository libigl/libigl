// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "fit_rotations.h"
#include "polar_svd3x3.h"
#include <igl/repmat.h>
#include <igl/verbose.h>
#include <igl/polar_dec.h>
#include <igl/polar_svd.h>
#include <iostream>

template <typename DerivedS, typename DerivedD>
IGL_INLINE void igl::fit_rotations(
  const Eigen::PlainObjectBase<DerivedS> & S,
        Eigen::PlainObjectBase<DerivedD> & R)
{
  using namespace std;
  const int dim = S.cols();
  const int nr = S.rows()/dim;
  assert(nr * dim == S.rows());

  // resize output
  R.resize(dim,dim*nr); // hopefully no op (should be already allocated)

  //std::cout<<"S=["<<std::endl<<S<<std::endl<<"];"<<std::endl;
  //MatrixXd si(dim,dim);
  Eigen::Matrix<typename DerivedS::Scalar,3,3> si;// = Eigen::Matrix3d::Identity();
  // loop over number of rotations we're computing
  for(int r = 0;r<nr;r++)
  {
    // build this covariance matrix
    for(int i = 0;i<dim;i++)
    {
      for(int j = 0;j<dim;j++)
      {
        si(i,j) = S(i*nr+r,j);
      }
    }
    Eigen::Matrix<typename DerivedD::Scalar,3,3> ri;
    Eigen::Matrix<typename DerivedD::Scalar,3,3> ti;
    //polar_dec(si,ri,ti);
    //polar_svd(si,ri,ti);
    polar_svd3x3(si, ri);
    assert(ri.determinant() >= 0);
#ifndef FIT_ROTATIONS_ALLOW_FLIPS
    // Check for reflection
    if(ri.determinant() < 0)
    {
      cerr<<"Error: Warning: flipping is wrong..."<<endl;
      assert(false && "This is wrong. Need to flip column in U and recompute R = U*V'");
      // flip sign of last row
      ri.row(2) *= -1;
    }
    assert(ri.determinant() >= 0);
#endif  
    // Not sure why polar_dec computes transpose...
    R.block(0,r*dim,dim,dim) = ri.block(0,0,dim,dim).transpose();
  }
}

template <typename DerivedS, typename DerivedD>
IGL_INLINE void igl::fit_rotations_planar(
  const Eigen::PlainObjectBase<DerivedS> & S,
        Eigen::PlainObjectBase<DerivedD> & R)
{ 
  using namespace std;
  const int dim = S.cols();
  const int nr = S.rows()/dim;
  assert(nr * dim == S.rows());

  // resize output
  R.resize(dim,dim*nr); // hopefully no op (should be already allocated)

  Eigen::Matrix<typename DerivedS::Scalar,2,2> si;
  // loop over number of rotations we're computing
  for(int r = 0;r<nr;r++)
  {
    // build this covariance matrix
    for(int i = 0;i<2;i++)
    {
      for(int j = 0;j<2;j++)
      {
        si(i,j) = S(i*nr+r,j);
      }
    }
    Eigen::Matrix<typename DerivedD::Scalar,2,2> ri;
    Eigen::Matrix<typename DerivedD::Scalar,2,2> ti;
    igl::polar_svd(si,ri,ti);
#ifndef FIT_ROTATIONS_ALLOW_FLIPS
    // Check for reflection
    if(ri.determinant() < 0)
    {
      cerr<<"Error: Warning: flipping is wrong..."<<endl;
      assert(false && "This is wrong. Need to flip column in U and recompute R = U*V'");
      // flip sign of last row
      ri.row(1) *= -1;
    }
    assert(ri.determinant() >= 0);
#endif  
    // Not sure why polar_dec computes transpose...
    R.block(0,r*dim,2,2) = ri.block(0,0,2,2).transpose();
    // Set remaining part to identity
    R(0,r*3+2) = 0;
    R(1,r*3+2) = 0;
    R(2,r*3+0) = 0;
    R(2,r*3+1) = 0;
    R(2,r*3+2) = 1;
  }
}


#ifdef __SSE__
IGL_INLINE void igl::fit_rotations_SSE(
  const Eigen::MatrixXf & S, 
  Eigen::MatrixXf & R)
{
  const int cStep = 4;

  assert(S.cols() == 3);
  const int dim = 3; //S.cols();
  const int nr = S.rows()/dim;  
  assert(nr * dim == S.rows());

  // resize output
  R.resize(dim,dim*nr); // hopefully no op (should be already allocated)

  Eigen::Matrix<float, 3*cStep, 3> siBig;
  // using SSE decompose cStep matrices at a time:
  int r = 0;
  for( ; r<nr; r+=cStep)
  {
    int numMats = cStep;
    if (r + cStep >= nr) numMats = nr - r;
    // build siBig:
    for (int k=0; k<numMats; k++)
    {
      for(int i = 0;i<dim;i++)
      {
        for(int j = 0;j<dim;j++)
        {
          siBig(i + 3*k, j) = S(i*nr + r + k, j);
        }
      }
    }
    Eigen::Matrix<float, 3*cStep, 3> ri;
    polar_svd3x3_sse(siBig, ri);    

    for (int k=0; k<cStep; k++)
      assert(ri.block(3*k, 0, 3, 3).determinant() >= 0);

    // Not sure why polar_dec computes transpose...
    for (int k=0; k<numMats; k++)
    {
      R.block(0, (r + k)*dim, dim, dim) = ri.block(3*k, 0, dim, dim).transpose();
    }    
  }
}

IGL_INLINE void igl::fit_rotations_SSE(
  const Eigen::MatrixXd & S,
  Eigen::MatrixXd & R)
{
  const Eigen::MatrixXf Sf = S.cast<float>();
  Eigen::MatrixXf Rf;
  fit_rotations_SSE(Sf,Rf);
  R = Rf.cast<double>();
}
#endif

#ifdef __AVX__
IGL_INLINE void igl::fit_rotations_AVX(
  const Eigen::MatrixXf & S,
  Eigen::MatrixXf & R)
{
  const int cStep = 8;

  assert(S.cols() == 3);
  const int dim = 3; //S.cols();
  const int nr = S.rows()/dim;  
  assert(nr * dim == S.rows());

  // resize output
  R.resize(dim,dim*nr); // hopefully no op (should be already allocated)

  Eigen::Matrix<float, 3*cStep, 3> siBig;
  // using SSE decompose cStep matrices at a time:
  int r = 0;
  for( ; r<nr; r+=cStep)
  {
    int numMats = cStep;
    if (r + cStep >= nr) numMats = nr - r;
    // build siBig:
    for (int k=0; k<numMats; k++)
    {
      for(int i = 0;i<dim;i++)
      {
        for(int j = 0;j<dim;j++)
        {
          siBig(i + 3*k, j) = S(i*nr + r + k, j);
        }
      }
    }
    Eigen::Matrix<float, 3*cStep, 3> ri;
    polar_svd3x3_avx(siBig, ri);    

    for (int k=0; k<cStep; k++)
      assert(ri.block(3*k, 0, 3, 3).determinant() >= 0);

    // Not sure why polar_dec computes transpose...
    for (int k=0; k<numMats; k++)
    {
      R.block(0, (r + k)*dim, dim, dim) = ri.block(3*k, 0, dim, dim).transpose();
    }    
  }
}
#endif

#ifndef IGL_HEADER_ONLY
// Explicit template instanciation
template void igl::fit_rotations<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::fit_rotations_planar<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::fit_rotations_planar<Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<float, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> >&);
#endif
