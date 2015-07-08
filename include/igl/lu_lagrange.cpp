// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "lu_lagrange.h"

// Cholesky LLT decomposition for symmetric positive definite
//#include <Eigen/SparseExtra>
// Bug in unsupported/Eigen/SparseExtra needs iostream first
#include <iostream>
#include <unsupported/Eigen/SparseExtra>
#include <cassert>
#include <cstdio>
#include "find.h"
#include "sparse.h"

template <typename T>
IGL_INLINE bool igl::lu_lagrange(
  const Eigen::SparseMatrix<T> & ATA,
  const Eigen::SparseMatrix<T> & C,
  Eigen::SparseMatrix<T> & L,
  Eigen::SparseMatrix<T> & U)
{
#if EIGEN_VERSION_AT_LEAST(3,0,92)
#if defined(_WIN32)
  #pragma message("lu_lagrange has not yet been implemented for your Eigen Version")
#else
  #warning lu_lagrange has not yet been implemented for your Eigen Version
#endif

  return false;
#else
  // number of unknowns
  int n = ATA.rows();
  // number of lagrange multipliers
  int m = C.cols();

  assert(ATA.cols() == n);
  if(m != 0)
  {
    assert(C.rows() == n);
    if(C.nonZeros() == 0)
    {
      // See note above about empty columns in C
      fprintf(stderr,"Error: lu_lagrange() C has columns but no entries\n");
      return false;
    }
  }

  // Check that each column of C has at least one entry
  std::vector<bool> has_entry; has_entry.resize(C.cols(),false);
  // Iterate over outside
  for(int k=0; k<C.outerSize(); ++k)
  {
    // Iterate over inside
    for(typename Eigen::SparseMatrix<T>::InnerIterator it (C,k); it; ++it)
    {
      has_entry[it.col()] = true;
    }
  }
  for(int i=0;i<(int)has_entry.size();i++)
  {
    if(!has_entry[i])
    {
      // See note above about empty columns in C
      fprintf(stderr,"Error: lu_lagrange() C(:,%d) has no entries\n",i);
      return false;
    }
  }



  // Cholesky factorization of ATA
  //// Eigen fails if you give a full view of the matrix like this:
  //Eigen::SparseLLT<SparseMatrix<T> > ATA_LLT(ATA);
  Eigen::SparseMatrix<T> ATA_LT = ATA.template triangularView<Eigen::Lower>();
  Eigen::SparseLLT<Eigen::SparseMatrix<T> > ATA_LLT(ATA_LT);

  Eigen::SparseMatrix<T> J = ATA_LLT.matrixL();

  //if(!ATA_LLT.succeeded())
  if(!((J*0).eval().nonZeros() == 0))
  {
    fprintf(stderr,"Error: lu_lagrange() failed to factor ATA\n");
    return false;
  }

  if(m == 0)
  {
    // If there are no constraints (C is empty) then LU decomposition is just L
    // and L' from cholesky decomposition
    L = J;
    U = J.transpose();
  }else
  {
    // Construct helper matrix M
    Eigen::SparseMatrix<T> M = C;
    J.template triangularView<Eigen::Lower>().solveInPlace(M);

    // Compute cholesky factorizaiton of M'*M
    Eigen::SparseMatrix<T> MTM = M.transpose() * M;

    Eigen::SparseLLT<Eigen::SparseMatrix<T> > MTM_LLT(MTM.template triangularView<Eigen::Lower>());

    Eigen::SparseMatrix<T> K = MTM_LLT.matrixL();

    //if(!MTM_LLT.succeeded())
    if(!((K*0).eval().nonZeros() == 0))
    {
      fprintf(stderr,"Error: lu_lagrange() failed to factor MTM\n");
      return false;
    }

    // assemble LU decomposition of Q
    Eigen::Matrix<int,Eigen::Dynamic,1> MI;
    Eigen::Matrix<int,Eigen::Dynamic,1> MJ;
    Eigen::Matrix<T,Eigen::Dynamic,1> MV;
    igl::find(M,MI,MJ,MV);

    Eigen::Matrix<int,Eigen::Dynamic,1> KI;
    Eigen::Matrix<int,Eigen::Dynamic,1> KJ;
    Eigen::Matrix<T,Eigen::Dynamic,1> KV;
    igl::find(K,KI,KJ,KV);

    Eigen::Matrix<int,Eigen::Dynamic,1> JI;
    Eigen::Matrix<int,Eigen::Dynamic,1> JJ;
    Eigen::Matrix<T,Eigen::Dynamic,1> JV;
    igl::find(J,JI,JJ,JV);

    int nnz = JV.size()  + MV.size() + KV.size();

    Eigen::Matrix<int,Eigen::Dynamic,1> UI(nnz);
    Eigen::Matrix<int,Eigen::Dynamic,1> UJ(nnz);
    Eigen::Matrix<T,Eigen::Dynamic,1> UV(nnz);
    UI << JJ,                        MI, (KJ.array() + n).matrix();
    UJ << JI, (MJ.array() + n).matrix(), (KI.array() + n).matrix(); 
    UV << JV,                        MV,                     KV*-1;
    igl::sparse(UI,UJ,UV,U);

    Eigen::Matrix<int,Eigen::Dynamic,1> LI(nnz);
    Eigen::Matrix<int,Eigen::Dynamic,1> LJ(nnz);
    Eigen::Matrix<T,Eigen::Dynamic,1> LV(nnz);
    LI << JI, (MJ.array() + n).matrix(), (KI.array() + n).matrix();
    LJ << JJ,                        MI, (KJ.array() + n).matrix(); 
    LV << JV,                        MV,                        KV;
    igl::sparse(LI,LJ,LV,L);
  }

  return true;
  #endif
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template bool igl::lu_lagrange<double>(Eigen::SparseMatrix<double, 0, int> const&, Eigen::SparseMatrix<double, 0, int> const&, Eigen::SparseMatrix<double, 0, int>&, Eigen::SparseMatrix<double, 0, int>&);
#endif
