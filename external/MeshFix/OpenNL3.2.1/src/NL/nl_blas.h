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

#ifndef __NL_BLAS__
#define __NL_BLAS__

#ifndef NL_FORTRAN_WRAP
#define NL_FORTRAN_WRAP(x) x##_
#endif


/***********************************************************************************/
/* C wrappers for BLAS routines */

/* x <- a*x */
void dscal( int n, double alpha, double *x, int incx ) ;

/* y <- x */
void dcopy( 
    int n, double *x, int incx, double *y, int incy 
) ;

/* y <- a*x+y */
void daxpy( 
    int n, double alpha, double *x, int incx, double *y,
    int incy 
) ;

/* returns x^T*y */
double ddot( 
    int n, double *x, int incx, double *y, int incy 
) ;

/** returns |x|_2 */
double dnrm2( int n, double *x, int incx ) ;

typedef enum { NoTranspose=0, Transpose=1, ConjugateTranspose=2 } MatrixTranspose ;
typedef enum { UpperTriangle=0, LowerTriangle=1 } MatrixTriangle ;
typedef enum { UnitTriangular=0, NotUnitTriangular=1 } MatrixUnitTriangular ;

/** x <- A^{-1}*x,  x <- A^{-T}*x */
void dtpsv( 
    MatrixTriangle uplo, MatrixTranspose trans,
    MatrixUnitTriangular diag, int n, double *AP,
    double *x, int incx 
) ;

/** y <- alpha*A*x + beta*y,  y <- alpha*A^T*x + beta*y,   A-(m,n) */
void dgemv( 
    MatrixTranspose trans, int m, int n, double alpha,
    double *A, int ldA, double *x, int incx,
    double beta, double *y, int incy 
) ;

#endif
