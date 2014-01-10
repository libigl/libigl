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

#include <NL/nl_matrix.h>

void nlRowColumnConstruct(NLRowColumn* c) {
    c->size     = 0 ;
    c->capacity = 0 ;
    c->coeff    = NULL ;
}

void nlRowColumnDestroy(NLRowColumn* c) {
    NL_DELETE_ARRAY(c->coeff) ;
#ifdef NL_PARANOID
    NL_CLEAR(NLRowColumn, c) ; 
#endif
}

void nlRowColumnGrow(NLRowColumn* c) {
    if(c->capacity != 0) {
        c->capacity = 2 * c->capacity ;
        c->coeff = NL_RENEW_ARRAY(NLCoeff, c->coeff, c->capacity) ;
    } else {
        c->capacity = 4 ;
        c->coeff = NL_NEW_ARRAY(NLCoeff, c->capacity) ;
    }
}

void nlRowColumnAdd(NLRowColumn* c, NLint index, NLdouble value) {
    NLuint i ;
    for(i=0; i<c->size; i++) {
        if(c->coeff[i].index == index) {
            c->coeff[i].value += value ;
            return ;
        }
    }
    if(c->size == c->capacity) {
        nlRowColumnGrow(c) ;
    }
    c->coeff[c->size].index = index ;
    c->coeff[c->size].value = value ;
    c->size++ ;
}

/* Does not check whether the index already exists */
void nlRowColumnAppend(NLRowColumn* c, NLint index, NLdouble value) {
    if(c->size == c->capacity) {
        nlRowColumnGrow(c) ;
    }
    c->coeff[c->size].index = index ;
    c->coeff[c->size].value = value ;
    c->size++ ;
}

void nlRowColumnZero(NLRowColumn* c) {
    c->size = 0 ;
}

void nlRowColumnClear(NLRowColumn* c) {
    c->size     = 0 ;
    c->capacity = 0 ;
    NL_DELETE_ARRAY(c->coeff) ;
}

static int nlCoeffCompare(const void* p1, const void* p2) {
    return (((NLCoeff*)(p2))->index > ((NLCoeff*)(p1))->index) ;
}

void nlRowColumnSort(NLRowColumn* c) {
    qsort(c->coeff, c->size, sizeof(NLCoeff), nlCoeffCompare) ;
}

/******************************************************************************/
/* SparseMatrix data structure */

void nlSparseMatrixConstruct(
    NLSparseMatrix* M, NLuint m, NLuint n, NLenum storage
) {
    NLuint i ;
    M->m = m ;
    M->n = n ;
    M->storage = storage ;
    if(storage & NL_MATRIX_STORE_ROWS) {
        M->row = NL_NEW_ARRAY(NLRowColumn, m) ;
        for(i=0; i<n; i++) {
            nlRowColumnConstruct(&(M->row[i])) ;
        }
    } else {
        M->row = NULL ;
    }

    if(storage & NL_MATRIX_STORE_COLUMNS) {
        M->column = NL_NEW_ARRAY(NLRowColumn, n) ;
        for(i=0; i<n; i++) {
            nlRowColumnConstruct(&(M->column[i])) ;
        }
    } else {
        M->column = NULL ;
    }

    M->diag_size = MIN(m,n) ;
    M->diag = NL_NEW_ARRAY(NLdouble, M->diag_size) ;
}

void nlSparseMatrixDestroy(NLSparseMatrix* M) {
    NLuint i ;
    NL_DELETE_ARRAY(M->diag) ;
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i=0; i<M->m; i++) {
            nlRowColumnDestroy(&(M->row[i])) ;
        }
        NL_DELETE_ARRAY(M->row) ;
    }
    if(M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i=0; i<M->n; i++) {
            nlRowColumnDestroy(&(M->column[i])) ;
        }
        NL_DELETE_ARRAY(M->column) ;
    }
#ifdef NL_PARANOID
    NL_CLEAR(NLSparseMatrix,M) ;
#endif
}

void nlSparseMatrixAdd(NLSparseMatrix* M, NLuint i, NLuint j, NLdouble value) {
    nl_parano_range_assert(i, 0, M->m - 1) ;
    nl_parano_range_assert(j, 0, M->n - 1) ;
    if((M->storage & NL_MATRIX_STORE_SYMMETRIC) && (j > i)) {
        return ;
    }
    if(i == j) {
        M->diag[i] += value ;
    }
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        nlRowColumnAdd(&(M->row[i]), j, value) ;
    }
    if(M->storage & NL_MATRIX_STORE_COLUMNS) {
        nlRowColumnAdd(&(M->column[j]), i, value) ;
    }
}

void nlSparseMatrixZero( NLSparseMatrix* M) {
    NLuint i ;
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i=0; i<M->m; i++) {
            nlRowColumnZero(&(M->row[i])) ;
        }
    }
    if(M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i=0; i<M->n; i++) {
            nlRowColumnZero(&(M->column[i])) ;
        }
    }
    NL_CLEAR_ARRAY(NLdouble, M->diag, M->diag_size) ;    
}

void nlSparseMatrixClear( NLSparseMatrix* M) {
    NLuint i ;
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i=0; i<M->m; i++) {
            nlRowColumnClear(&(M->row[i])) ;
        }
    }
    if(M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i=0; i<M->n; i++) {
            nlRowColumnClear(&(M->column[i])) ;
        }
    }
    NL_CLEAR_ARRAY(NLdouble, M->diag, M->diag_size) ;    
}

/* Returns the number of non-zero coefficients */
NLuint nlSparseMatrixNNZ( NLSparseMatrix* M) {
    NLuint nnz = 0 ;
    NLuint i ;
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i = 0; i<M->m; i++) {
            nnz += M->row[i].size ;
        }
    } else if (M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i = 0; i<M->n; i++) {
            nnz += M->column[i].size ;
        }
    } else {
        nl_assert_not_reached ;
    }
    return nnz ;
}

void nlSparseMatrixSort( NLSparseMatrix* M) {
    NLuint i ;
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i = 0; i<M->m; i++) {
            nlRowColumnSort(&(M->row[i])) ;
        }
    } 
    if (M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i = 0; i<M->n; i++) {
            nlRowColumnSort(&(M->row[i])) ;
        }
    } 
}

/************************************************************************************/
/* SparseMatrix x Vector routines, internal helper routines */

static void nlSparseMatrix_mult_rows_symmetric(
        NLSparseMatrix* A,
        NLdouble* x,
        NLdouble* y) {
    NLuint m = A->m ;
    NLuint i,ij ;
    NLRowColumn* Ri = NULL ;
    NLCoeff* c = NULL ;
    for(i=0; i<m; i++) {
        y[i] = 0 ;
        Ri = &(A->row[i]) ;
        for(ij=0; ij<Ri->size; ij++) {
            c = &(Ri->coeff[ij]) ;
            y[i] += c->value * x[c->index] ;
            if(i != c->index) {
                y[c->index] += c->value * x[i] ;
            }
        }
    }
}

static void nlSparseMatrix_mult_rows(
        NLSparseMatrix* A,
        NLdouble* x,
        NLdouble* y) {
    NLuint m = A->m ;
    NLuint i,ij ;
    NLRowColumn* Ri = NULL ;
    NLCoeff* c = NULL ;
    for(i=0; i<m; i++) {
        y[i] = 0 ;
        Ri = &(A->row[i]) ;
        for(ij=0; ij<Ri->size; ij++) {
            c = &(Ri->coeff[ij]) ;
            y[i] += c->value * x[c->index] ;
        }
    }
}

static void nlSparseMatrix_mult_cols_symmetric(
        NLSparseMatrix* A,
        NLdouble* x,
        NLdouble* y) {
    NLuint n = A->n ;
    NLuint j,ii ;
    NLRowColumn* Cj = NULL ;
    NLCoeff* c = NULL ;
    for(j=0; j<n; j++) {
        y[j] = 0 ;
        Cj = &(A->column[j]) ;
        for(ii=0; ii<Cj->size; ii++) {
            c = &(Cj->coeff[ii]) ;
            y[c->index] += c->value * x[j] ;
            if(j != c->index) {
                y[j] += c->value * x[c->index] ;
            }
        }
    }
}

static void nlSparseMatrix_mult_cols(
        NLSparseMatrix* A,
        NLdouble* x,
        NLdouble* y) {
    NLuint n = A->n ;
    NLuint j,ii ; 
    NLRowColumn* Cj = NULL ;
    NLCoeff* c = NULL ;
    NL_CLEAR_ARRAY(NLdouble, y, A->m) ;
    for(j=0; j<n; j++) {
        Cj = &(A->column[j]) ;
        for(ii=0; ii<Cj->size; ii++) {
            c = &(Cj->coeff[ii]) ;
            y[c->index] += c->value * x[j] ;
        }
    }
}

/************************************************************************************/
/* SparseMatrix x Vector routines, main driver routine */

void nlSparseMatrixMult(NLSparseMatrix* A, NLdouble* x, NLdouble* y) {
    if(A->storage & NL_MATRIX_STORE_ROWS) {
        if(A->storage & NL_MATRIX_STORE_SYMMETRIC) {
            nlSparseMatrix_mult_rows_symmetric(A, x, y) ;
        } else {
            nlSparseMatrix_mult_rows(A, x, y) ;
        }
    } else {
        if(A->storage & NL_MATRIX_STORE_SYMMETRIC) {
            nlSparseMatrix_mult_cols_symmetric(A, x, y) ;
        } else {
            nlSparseMatrix_mult_cols(A, x, y) ;
        }
    }
}

