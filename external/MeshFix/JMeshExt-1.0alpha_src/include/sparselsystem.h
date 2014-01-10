/****************************************************************************
* JMeshExt                                                                  *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
*                                                                           *
* Copyright(C) 2006: IMATI-GE / CNR                                         *
*                                                                           *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifndef SPARSELSYSTEM_H
#define SPARSELSYSTEM_H

#include "matrix.h"

//////////////////////////////////////////////////////////////////////////
//
// Sparse linear system
//
//////////////////////////////////////////////////////////////////////////

//! Sparse linear system Ax = B

//! Matrix A is a square matrix of size 'system_size'.
//! A is initially 0. Non-zero coefficients can be summed through
//! the method 'sumCoefficient(coeff_value, row, column)'.
//! Matrix B is a generic matrix having 'system_size' rows and 'kterm_size' columns.
//! B is initially 0. B[i][j] can be set by summing
//! values through sumKnownTerm(value, row, column).
//! A solution of the system can be obtained using one column of B trhough
//! solve(result, which_column_of_B_to_use).
//! In most geometric algorithms 'kterm_size' is 3, and the searched solution
//! is the position of some points. In this case, the value of the coordinates
//! may be retrieved by solving the system three times, once to retrieve the 'x'
//! coordinates, once for the 'y's and once for the 'z's, using the first, the
//! second and the third column of B respectively.

class sparseSystem
{
 protected:

 class coeffIndexPair
 {
  public:
  int index;
  double coeff;

  coeffIndexPair(int a, double b) {index=a; coeff=b;}
 };

 class sparseSystemRow
 {
  public:

  List cips;
  sparseSystemRow() {}
  ~sparseSystemRow() {cips.freeNodes();}

  void addCoefficient(int, double);
  void print(FILE *, int);

  static int rowcompare(const void *, const void *);
 };

 protected:

 int num_equations;	//!< Number of equations (rows)
 int num_variables;	//!< Number of variables (columns)
 int kterm_size;	//!< Nr. of columns of the known term
 sparseSystemRow *rows;	//!< Rows of the system
 double **known_term;	//!< Actual coefficients of the known term

 public:

 sparseSystem(int s, int k, int n = 0); //!< Constructs an s x n system having a k column-wide known term
 ~sparseSystem();		//!< Destructor

 void sumCoefficient(double v, int i, int j) {rows[i].addCoefficient(j, v);} //!< Sums 'v' to A[i][j]
 void setKnownTerm(double v, int i, int j) {known_term[j][i] =v;} //!< Sets j'th component of B[i] to 'v' (B[i][j]=v)
 void sumKnownTerm(double v, int i, int j) {known_term[j][i]+=v;} //!< Sums j'th component of B[i] to 'v' (B[i][j]=v)
 bool solve(double *solution, int j); //!< Solves the system for j'th component of B (j starts from 0). False if fails.
 void print(FILE * =stdout);
};

// Sparse linear system with size-3 veriables
class sparse3System : public sparseSystem
{
 bool *locks;

 public:

 sparse3System(int s) : sparseSystem(s, 3)
  {locks = new bool[num_variables]; for (int i=0; i<num_variables; i++) locks[i]=false;}
 ~sparse3System() {delete locks;}

 // The following calls 'solve' three times and fills the vertex array with x, y, z results
 void solve(double *);
 void sumKnownTerm(const double *v, int i);
 void lock(int i) {locks[i]=true;}
 void unlock(int i) {locks[i]=false;}
};

// Underdetermined sparse system with size-1 veriables to be solved in the least-squares sense
class leastSquaresSystem : public sparseSystem
{
 bool *locks;

 public:

 leastSquaresSystem(int s, int n = 0);
 ~leastSquaresSystem() {delete locks;}

 void solve(double *);
 void setKnownTerm(const double v, int i) {sparseSystem::setKnownTerm(v, i, 0);}
 void sumKnownTerm(const double v, int i) {sparseSystem::sumKnownTerm(v, i, 0);}
 void lock(int i) {locks[i]=true;}
 void unlock(int i) {locks[i]=false;}
};

#endif // SPARSELSYSTEM_H
