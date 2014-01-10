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

#ifndef __NL_CONTEXT__
#define __NL_CONTEXT__

#include <NL/nl_private.h>
#include <NL/nl_matrix.h>

/******************************************************************************/
/* NLContext data structure */

typedef void(*NLMatrixFunc)(double* x, double* y) ;
typedef NLboolean(*NLSolverFunc)() ;

typedef struct {
    NLdouble  value ;
    NLboolean locked ;
    NLuint    index ;
} NLVariable ;

#define NL_STATE_INITIAL                0
#define NL_STATE_SYSTEM                 1
#define NL_STATE_MATRIX                 2
#define NL_STATE_ROW                    3
#define NL_STATE_MATRIX_CONSTRUCTED     4
#define NL_STATE_SYSTEM_CONSTRUCTED     5
#define NL_STATE_SOLVED                 6

typedef struct {
    NLenum           state ;
    NLVariable*      variable ;
    NLuint           n ;
    NLSparseMatrix   M ;
    NLRowColumn      af ;
    NLRowColumn      al ;
    NLRowColumn      xl ;
    NLdouble*        x ;
    NLdouble*        b ;
    NLdouble         right_hand_side ;
    NLdouble         row_scaling ;
    NLenum           solver ;
    NLenum           preconditioner ;
    NLuint           nb_variables ;
    NLuint           current_row ;
    NLboolean        least_squares ;
    NLboolean        symmetric ;
    NLuint           max_iterations ;
    NLuint           inner_iterations ;
    NLdouble         threshold ;
    NLdouble         omega ;
    NLboolean        normalize_rows ;
    NLboolean        alloc_M ;
    NLboolean        alloc_af ;
    NLboolean        alloc_al ;
    NLboolean        alloc_xl ;
    NLboolean        alloc_variable ;
    NLboolean        alloc_x ;
    NLboolean        alloc_b ;
    NLuint           used_iterations ;
    NLdouble         error ;
    NLdouble         elapsed_time ;
    NLMatrixFunc     matrix_vector_prod ;
    NLMatrixFunc     precond_vector_prod ;
    NLSolverFunc     solver_func ;
} NLContextStruct ;

extern NLContextStruct* nlCurrentContext ;

void nlCheckState(NLenum state) ;
void nlTransition(NLenum from_state, NLenum to_state) ;

void nlMatrixVectorProd_default(NLdouble* x, NLdouble* y) ;
NLboolean nlDefaultSolver() ;

#endif
