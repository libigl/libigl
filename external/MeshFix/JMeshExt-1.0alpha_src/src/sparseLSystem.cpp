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

#include "sparseLSystem.h"
#include "nl.h"

//////////////////////////////////////////////////////////////////////////
//
// Sparse linear system
//
//////////////////////////////////////////////////////////////////////////

void sparseSystem::sparseSystemRow::addCoefficient(int i, double c)
{
 Node *n;
 coeffIndexPair *f;
 for (n=cips.head(); n!=NULL; n=n->next())
 {
  f = (coeffIndexPair *)n->data;
  if (f->index == i) {f->coeff += c; return;}
 }
 cips.appendTail(new coeffIndexPair(i,c));
}

int sparseSystem::sparseSystemRow::rowcompare(const void *a, const void *b)
{
 coeffIndexPair *p1 = (coeffIndexPair *)a;
 coeffIndexPair *p2 = (coeffIndexPair *)b;

 if (p1->index > p2->index) return 1;
 if (p1->index < p2->index) return -1;
 return 0;
}

void sparseSystem::sparseSystemRow::print(FILE *fp, int sz)
{
 List stc(cips);
 stc.sort(&sparseSystem::sparseSystemRow::rowcompare);

 coeffIndexPair *f;
 int j=0;
 for (Node *n=stc.head(); n!=NULL; n=n->next())
 {
  f = (coeffIndexPair *)n->data;
  while (f->index > j) {fprintf(fp,"0.000 "); j++;}
  fprintf(fp, "%.3f ",f->coeff); j++;
 }
 for (; j<sz; j++) fprintf(fp,"0.000 ");
}

sparseSystem::sparseSystem(int s, int k, int n)
{
 num_equations = (n==0)?(s):(n);
 num_variables = s;
 kterm_size = k;
 rows = new sparseSystemRow[num_equations];
 known_term = new double *[kterm_size];
 int i,j;
 for (i=0; i<kterm_size; i++)
 {
  known_term[i] = new double[num_equations];
  for (j=0; j<num_equations; j++) known_term[i][j]=0.0;
 }
}

sparseSystem::~sparseSystem()
{
 delete [] rows;
 delete known_term;
}


// Solves the system for j'th component of B
bool sparseSystem::solve(double *x, int j)
{
 int i;
 Node *n;
 coeffIndexPair *f;

 nlNewContext();
 nlSolverParameteri(NL_SOLVER, NL_PERM_SUPERLU_EXT);
 nlSolverParameteri(NL_NB_VARIABLES, num_variables);
 nlBegin(NL_SYSTEM);

 for (i=0; i<num_variables; i++) nlSetVariable(i, x[i]);

 nlBegin(NL_MATRIX);

 for (i=0; i<num_equations; i++)
 {
  nlRowParameterd(NL_RIGHT_HAND_SIDE, known_term[j][i]);
  nlBegin(NL_ROW);
  for (n=rows[i].cips.head(); n!=NULL; n=n->next())
  {
   f = (coeffIndexPair *)n->data;
   nlCoefficient(f->index, f->coeff);
  }  
  nlEnd(NL_ROW);
 }

 nlEnd(NL_MATRIX);
 nlEnd(NL_SYSTEM);
 bool success = (bool)nlSolve();

 if (success) for (i=0; i<num_variables; i++) x[i] = nlGetVariable(i);

 nlDeleteContext(nlGetCurrent());

 return success;
}

void sparseSystem::print(FILE *fp)
{
 int j;

 for (int i=0; i<num_equations; i++)
 {
  rows[i].print(fp, num_variables);
  fprintf(fp,": ");
  for (j=0; j<kterm_size; j++) fprintf(fp,"%.3f ",known_term[j][i]);
  fprintf(fp,"\n");
 }
}


void sparse3System::solve(double *vs)
{
 int i, j;
 Node *n;
 coeffIndexPair *f;

 for (j=0; j<3; j++)
 {
  nlNewContext();
  nlSolverParameteri(NL_SOLVER, NL_PERM_SUPERLU_EXT);
  nlSolverParameteri(NL_NB_VARIABLES, num_variables);
  nlBegin(NL_SYSTEM);

  for (i=0; i<num_variables; i++)
  {
   nlSetVariable(i, vs[i*3 + j]);
   if (locks[i]) nlLockVariable(i);
  }

  nlBegin(NL_MATRIX);

  for (i=0; i<num_equations; i++)
  {
   nlRowParameterd(NL_RIGHT_HAND_SIDE, known_term[j][i]);
   nlBegin(NL_ROW);
   for (n=rows[i].cips.head(); n!=NULL; n=n->next())
   {
    f = (coeffIndexPair *)n->data;
    nlCoefficient(f->index, f->coeff);
   }  
   nlEnd(NL_ROW);
  }

  nlEnd(NL_MATRIX);
  nlEnd(NL_SYSTEM);
  nlSolve();

  for (i=0; i<num_variables; i++) vs[i*3 + j] = nlGetVariable(i);

  nlDeleteContext(nlGetCurrent());
 }
}


void sparse3System::sumKnownTerm(const double *v, int i)
{
 sparseSystem::sumKnownTerm(v[0], i, 0);
 sparseSystem::sumKnownTerm(v[1], i, 1);
 sparseSystem::sumKnownTerm(v[2], i, 2);
}


leastSquaresSystem::leastSquaresSystem(int s, int n) : sparseSystem(s, 1, n)
{
 locks = new bool[num_variables];
 for (int i=0; i<num_variables; i++) locks[i]=false;
}

void leastSquaresSystem::solve(double *vs)
{
 int i;
 Node *n;
 coeffIndexPair *f;

 nlNewContext();
 nlSolverParameteri(NL_SOLVER, NL_PERM_SUPERLU_EXT);
 nlSolverParameteri(NL_NB_VARIABLES, num_variables);
 nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE) ;
 nlSolverParameteri(NL_MAX_ITERATIONS, 1000) ;
 nlSolverParameterd(NL_THRESHOLD, 1e-10) ;

 nlBegin(NL_SYSTEM);

 for (i=0; i<num_variables; i++)
 {
  nlSetVariable(i, vs[i]);
  if (locks[i]) nlLockVariable(i);
 }

 nlBegin(NL_MATRIX);

 for (i=0; i<num_equations; i++)
 {
  nlRowParameterd(NL_RIGHT_HAND_SIDE, known_term[0][i]);
  nlBegin(NL_ROW);
  for (n=rows[i].cips.head(); n!=NULL; n=n->next())
  {
   f = (coeffIndexPair *)n->data;
   nlCoefficient(f->index, f->coeff);
  }  
  nlEnd(NL_ROW);
 }

 nlEnd(NL_MATRIX);
 nlEnd(NL_SYSTEM);
 nlSolve();

 for (i=0; i<num_variables; i++) vs[i] = nlGetVariable(i);

 nlDeleteContext(nlGetCurrent());
}
