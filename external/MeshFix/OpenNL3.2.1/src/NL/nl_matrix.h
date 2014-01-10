/*
 *  Copyright (c) 2004-2010, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     levy@loria.fr
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#include <NL/nl_private.h>

#ifndef __NL_MATRIX__
#define __NL_MATRIX__

#ifdef __cplusplus
extern "C" {
#endif

/************************************************************************************/
/* Dynamic arrays for sparse row/columns */

typedef struct  {
    NLuint index ;
    NLdouble value ;
} NLCoeff ;

typedef struct {
    NLuint size ;
    NLuint capacity ;
    NLCoeff* coeff ;
} NLRowColumn ;

void nlRowColumnConstruct(NLRowColumn* c) ;
void nlRowColumnDestroy(NLRowColumn* c) ;
void nlRowColumnGrow(NLRowColumn* c) ;
void nlRowColumnAdd(NLRowColumn* c, NLint index, NLdouble value) ;
void nlRowColumnAppend(NLRowColumn* c, NLint index, NLdouble value) ;
void nlRowColumnZero(NLRowColumn* c) ;
void nlRowColumnClear(NLRowColumn* c) ;
void nlRowColumnSort(NLRowColumn* c) ;

/************************************************************************************/
/* SparseMatrix data structure */

#define NL_MATRIX_STORE_ROWS      1
#define NL_MATRIX_STORE_COLUMNS   2
#define NL_MATRIX_STORE_SYMMETRIC 4

typedef struct {
    NLuint m ;
    NLuint n ;
    NLuint diag_size ;
    NLenum storage ;
    NLRowColumn* row ;
    NLRowColumn* column ;
    NLdouble*      diag ;
} NLSparseMatrix ;


void nlSparseMatrixConstruct(
    NLSparseMatrix* M, NLuint m, NLuint n, NLenum storage
) ;

void nlSparseMatrixDestroy(NLSparseMatrix* M) ;

void nlSparseMatrixAdd(
    NLSparseMatrix* M, NLuint i, NLuint j, NLdouble value
) ;

void nlSparseMatrixZero( NLSparseMatrix* M) ;
void nlSparseMatrixClear( NLSparseMatrix* M) ;
NLuint nlSparseMatrixNNZ( NLSparseMatrix* M) ;
void nlSparseMatrixSort( NLSparseMatrix* M) ;

/************************************************************************************/
/* SparseMatrix x Vector routine */

void nlSparseMatrixMult(NLSparseMatrix* A, NLdouble* x, NLdouble* y) ;

#ifdef __cplusplus
}
#endif

#endif
