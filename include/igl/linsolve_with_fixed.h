// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_LINSOLVE_WITH_FIXED_H
#define IGL_LINSOLVE_WITH_FIXED_H
#include "igl_inline.h"

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseExtra>
// Bug in unsupported/Eigen/SparseExtra needs iostream first
#include <iostream>
#include <unsupported/Eigen/SparseExtra>

namespace igl
{
  // LINSOLVE_WITH_FIXED Solves a sparse linear system Ax=b with
  // fixed known values in x

  // Templates:
  //   T  should be a eigen matrix primitive type like float or double
  // Inputs:
  //   A  m by n sparse matrix
  //   b  n by k dense matrix, k is bigger than 1 for multiple right handsides
  //   known list of indices to known rows in x (must be sorted in ascending order)
  //   Y  list of fixed values corresponding to known rows in Z
  // Outputs:
  //   x  n by k dense matrix containing the solution

  template<typename T, typename Derived>
  inline void SolveLinearSystemWithBoundaryConditions(const SparseMatrix<T>& A, const MatrixBase<Derived>& b, const vector<int>& known, MatrixBase<Derived>& x)
  {
  const int rows = x.rows();
  const int cols = x.cols();
  SparseMatrix<T> AP(rows,rows);

  const int knownCount = known.size();
  const int unknownCount = rows - knownCount;
  int knownIndex = 0;
  int unknownIndex = 0;

  if(knownCount >= rows)
  {
  std::cerr << "No unknowns to solve for!\n";
  return;
  }

  // Permutate to sort by known boundary conditions -> [A1,A2] * [x1;x2]
  VectorXi transpositions;
  transpositions.resize(rows);
  for(int n=0;n<rows;n++)
  {
    if(knownIndex < known.size() && n == known[knownIndex])
    {
      transpositions[unknownCount + knownIndex] = n;
      knownIndex++;
    }
    else
    {
      transpositions[unknownIndex] = n;
      unknownIndex++;
    }
  }
  PermutationMatrix<Dynamic,Dynamic,int> P(transpositions);
  PermutationMatrix<Dynamic,Dynamic,int> PInv = P.inverse();
  A.twistedBy(PInv).evalTo(AP);

  // Split in kown and unknown parts A1,A2,x2,b1
  SparseMatrix<T> A1, A2;
  internal::plain_matrix_type<Derived>::type x1, x2, b1;
  x2 = (PInv*x).block(unknownCount,0,knownCount,cols);
  b1 = (PInv*b).block(0,0,unknownCount,cols);

  SparseMatrix<T> A1Temp(rows,unknownCount);
  SparseMatrix<T> A2Temp(rows,knownCount);
  A1Temp = AP.innerVectors(0,unknownCount).transpose();
  A2Temp = AP.innerVectors(unknownCount,knownCount).transpose();
  A1 = A1Temp.innerVectors(0,unknownCount).transpose();
  A2 = A2Temp.innerVectors(0,unknownCount).transpose();

  // Solve A1*x1 = b - A2*x2
  SparseLU<SparseMatrix<T>> solver;
  solver.compute(A1);
  if (solver.info() != Eigen::Success)
  {
  std::cerr << "Matrix can not be decomposed.\n";
  return;
  }

  internal::plain_matrix_type<Derived>::type t = b1-A2*x2;
  x1 = solver.solve(t);
  if (solver.info() != Eigen::Success)
  {
  std::cerr << "Solving failed.\n";
  return;
  }

  // Compose resulting x
  x.block(0,0,unknownCount,cols) = x1;
  x.block(unknownCount,0,knownCount,cols) = x2;
  x = P * x;
  };
}