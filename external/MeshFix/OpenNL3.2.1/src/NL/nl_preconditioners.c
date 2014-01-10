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

#include <NL/nl_preconditioners.h>
#include <NL/nl_blas.h>
#include <NL/nl_matrix.h>
#include <NL/nl_context.h>

/******************************************************************************/
/* preconditioners */

/* Utilities for preconditioners */

void nlMultDiagonal(NLdouble* xy, NLdouble omega) {
    NLuint N = nlCurrentContext->n ;
    NLuint i ;
    NLdouble* diag = nlCurrentContext->M.diag ;
    for(i=0; i<N; i++) {
        xy[i] *= (diag[i] / omega) ;
    }
}

void nlMultDiagonalInverse(NLdouble* xy, NLdouble omega) {
    NLuint N = nlCurrentContext->n ;
    NLuint i ;
    NLdouble* diag = nlCurrentContext->M.diag ;
    for(i=0; i<N; i++) {
	    xy[i] *= ((diag[i] != 0) ? (omega / diag[i]) : omega) ;
    }
}

void nlMultLowerInverse(NLdouble* x, NLdouble* y, double omega) {
    NLSparseMatrix* A = &(nlCurrentContext->M) ;
    NLuint n       = A->n ;
    NLdouble* diag = A->diag ;
    NLuint i ;
    NLuint ij ;
    NLRowColumn*  Ri = NULL ;
    NLCoeff* c = NULL ;
    NLdouble S ;

    nl_assert(A->storage & NL_MATRIX_STORE_SYMMETRIC) ;
    nl_assert(A->storage & NL_MATRIX_STORE_ROWS) ;

    for(i=0; i<n; i++) {
        S = 0 ;
        Ri = &(A->row[i]) ;
        for(ij=0; ij < Ri->size; ij++) {
            c = &(Ri->coeff[ij]) ;
            nl_parano_assert(c->index <= i) ; 
            if(c->index != i) {
                S += c->value * y[c->index] ; 
            }
        }
        y[i] = (x[i] - S) * omega / diag[i] ;
    }
}

void nlMultUpperInverse(NLdouble* x, NLdouble* y, NLdouble omega) {
    NLSparseMatrix* A = &(nlCurrentContext->M) ;
    NLuint n       = A->n ;
    NLdouble* diag = A->diag ;
    NLint i ;
    NLuint ij ;
    NLRowColumn*  Ci = NULL ;
    NLCoeff* c = NULL ;
    NLdouble S ;

    nl_assert(A->storage & NL_MATRIX_STORE_SYMMETRIC) ;
    nl_assert(A->storage & NL_MATRIX_STORE_COLUMNS) ;

    for(i=n-1; i>=0; i--) {
        S = 0 ;
        Ci = &(A->column[i]) ;
        for(ij=0; ij < Ci->size; ij++) {
            c = &(Ci->coeff[ij]) ;
            nl_parano_assert(c->index >= i) ; 
            if(c->index != i) {
                S += c->value * y[c->index] ; 
            }
        }
        y[i] = (x[i] - S) * omega / diag[i] ;
    }
}


void nlPreconditioner_Jacobi(NLdouble* x, NLdouble* y) {
    NLuint N = nlCurrentContext->n ;
    dcopy(N, x, 1, y, 1) ;
    nlMultDiagonalInverse(y, 1.0) ;
}

void nlPreconditioner_SSOR(NLdouble* x, NLdouble* y) {
    NLdouble omega = nlCurrentContext->omega ;
    static double* work = NULL ;
    static int work_size = 0 ;
    NLuint n = nlCurrentContext->n ;
    if(n != work_size) {
        work = NL_RENEW_ARRAY(NLdouble, work, n) ;
        work_size = n ;
    }
    
    nlMultLowerInverse(x, work, omega) ;
    nlMultDiagonal(work, omega) ;
    nlMultUpperInverse(work, y, omega) ;

    dscal(n, 2.0 - omega, y, 1) ;
}

