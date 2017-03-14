// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "eigs.h"

#include "cotmatrix.h"
#include "sort.h"
#include "slice.h"
#include "massmatrix.h"
#include <iostream>

template <
  typename Atype,
  typename Btype,
  typename DerivedU,
  typename DerivedS>
IGL_INLINE bool igl::eigs(
  const Eigen::SparseMatrix<Atype> & A,
  const Eigen::SparseMatrix<Btype> & iB,
  const EigsType type,
  Eigen::PlainObjectBase<DerivedU> & sU,
  Eigen::PlainObjectBase<DerivedS> & sS,
  const size_t k,
  unsigned int max_iter)
{
  using namespace Eigen;
  using namespace std;
  const size_t n = A.rows();
  assert(A.cols() == n && "A should be square.");
  assert(iB.rows() == n && "B should be match A's dims.");
  assert(iB.cols() == n && "B should be square.");
  assert(type == EIGS_TYPE_SM && "Only low frequencies are supported");
  DerivedU U(n, 1);
  DerivedS S(1, 1);
  typedef Atype Scalar;
  typedef Eigen::Matrix<typename DerivedU::Scalar,DerivedU::RowsAtCompileTime,1> VectorXS;
  // Rescale B for better numerics
  const Scalar rescale = std::abs(iB.diagonal().maxCoeff());
  const Eigen::SparseMatrix<Btype> B = iB/rescale;

  Scalar tol = 1e-4;
  Scalar conv = 1e-14;
  int i = 0;
  while(true)
  {
    // Random initial guess
    VectorXS y = VectorXS::Random(n,1);
    Scalar eff_sigma = 0;
    if(i>0)
    {
      eff_sigma = 1e-8+std::abs(S(i-1));
    }
    // whether to use rayleigh quotient method
    bool ray = false;
    Scalar err = std::numeric_limits<Scalar>::infinity();
    int iter;
    Scalar sigma = std::numeric_limits<Scalar>::infinity();
    VectorXS x;
    for(iter = 0;iter<max_iter;iter++)
    {
      if(i>0 && !ray)
      {
        // project-out existing modes
        for(int j = 0;j<i;j++)
        {
          const VectorXS u = U.col(j);
          y = (y - u*u.dot(B*y)/u.dot(B * u)).eval();
        }
      }
      // normalize
      x = y/sqrt(y.dot(B*y));

      // current guess at eigen value
      sigma = x.dot(A*x)/x.dot(B*x);
      //x *= sigma>0?1.:-1.;

      Scalar err_prev = err;
      err = (A*x-sigma*B*x).array().abs().maxCoeff();
      if(err<conv)
      {
        break;
      }
      if(ray || err<tol)
      {
        eff_sigma = sigma;
        ray = true;
      }

      Scalar tikhonov = std::abs(eff_sigma)<1e-12?1e-10:0;
      switch(type)
      {
        default:
          assert(false && "Not supported");
          break;
        case EIGS_TYPE_SM:
        {
          SimplicialLDLT<SparseMatrix<Scalar> > solver;
          const SparseMatrix<Scalar> C = A-eff_sigma*B+tikhonov*B;
          //mw.save(C,"C");
          //mw.save(eff_sigma,"eff_sigma");
          //mw.save(tikhonov,"tikhonov");
          solver.compute(C);
          switch(solver.info())
          {
            case Eigen::Success:
              break;
            case Eigen::NumericalIssue:
              cerr<<"Error: Numerical issue."<<endl;
              return false;
            default:
              cerr<<"Error: Other."<<endl;
              return false;
          }
          const VectorXS rhs = B*x;
          y = solver.solve(rhs);
          //mw.save(rhs,"rhs");
          //mw.save(y,"y");
          //mw.save(x,"x");
          //mw.write("eigs.mat");
          //if(i == 1)
          //return false;
          break;
        }
      }
    }
    if(iter == max_iter)
    {
	  cout << "break by convergence failure on iteration [ " << i << " ]" << std::endl;
	  break;
    }
    if(i==0 || (S.head(i).array()-sigma).abs().maxCoeff()>1e-14)
    {
      U.col(i) = x;
      S(i) = sigma;
      i++;

	  if (k != 0 && i >= k) {
		std::cout << "found all eVec/eVals defined by k = " << k << std::endl;
		break;
	  }

	  U.conservativeResize(n, i+1);
	  S.conservativeResize(i+1, 1);
    }else
    {
      // restart with new random guess.
      cout<<"RESTART!"<<endl;
    }
  }
  // finally sort
  VectorXi I;
  igl::sort(S,1,false,sS,I);
  sU = igl::slice(U,I,2);
  sS /= rescale;
  sU /= sqrt(rescale);
  return true;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template bool igl::eigs<double, double, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::SparseMatrix<double, 0, int> const&, Eigen::SparseMatrix<double, 0, int> const&, igl::EigsType, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, unsigned long, unsigned int);
#ifdef WIN32
template bool igl::eigs<double, double, Eigen::Matrix<double,-1,-1,0,-1,-1>, Eigen::Matrix<double,-1,1,0,-1,1> >(Eigen::SparseMatrix<double,0,int> const &,Eigen::SparseMatrix<double,0,int> const &, igl::EigsType, Eigen::PlainObjectBase< Eigen::Matrix<double,-1,-1,0,-1,-1> > &, Eigen::PlainObjectBase<Eigen::Matrix<double,-1,1,0,-1,1> > &, unsigned long long, unsigned int);
#endif
#endif
