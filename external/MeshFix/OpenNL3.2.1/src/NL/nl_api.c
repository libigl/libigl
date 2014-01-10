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
#include <NL/nl_matrix.h>
#include <NL/nl_context.h>
#include <NL/nl_iterative_solvers.h>
#include <NL/nl_preconditioners.h>
#include <NL/nl_cnc_gpu_cuda.h>
#include <NL/nl_superlu.h>

NLboolean nlInitExtension(const char* extension) {
#ifdef NL_USE_SUPERLU
    if(!strcmp(extension, "SUPERLU")) {
        return NL_TRUE ;
    }
#endif
#ifdef NL_USE_CNC
    if(!strcmp(extension, "CNC")) {
        return NL_TRUE ;
    }
#endif
    return NL_FALSE ;
}



/************************************************************************************/
/* Get/Set parameters */

void nlSolverParameterd(NLenum pname, NLdouble param) {
    nlCheckState(NL_STATE_INITIAL) ;
    switch(pname) {
    case NL_SOLVER: {
        nlCurrentContext->solver = (NLenum)param ;
    } break ;
    case NL_NB_VARIABLES: {
        nl_assert(param > 0) ;
        nlCurrentContext->nb_variables = (NLuint)param ;
    } break ;
    case NL_LEAST_SQUARES: {
        nlCurrentContext->least_squares = (NLboolean)param ;
    } break ;
    case NL_MAX_ITERATIONS: {
        nl_assert(param > 0) ;
        nlCurrentContext->max_iterations = (NLuint)param ;
    } break ;
    case NL_THRESHOLD: {
        nl_assert(param >= 0) ;
        nlCurrentContext->threshold = (NLdouble)param ;
    } break ;
    case NL_OMEGA: {
        nl_range_assert(param,1.0,2.0) ;
        nlCurrentContext->omega = (NLdouble)param ;
    } break ;
    case NL_SYMMETRIC: {
        nlCurrentContext->symmetric = (NLboolean)param ;        
    }
    case NL_INNER_ITERATIONS: {
        nl_assert(param > 0) ;
        nlCurrentContext->inner_iterations = (NLuint)param ;
    } break ;
    case NL_PRECONDITIONER: {
        nlCurrentContext->preconditioner = (NLuint)param ;        
    } break ;
    default: {
        nl_assert_not_reached ;
    } break ;
    }
}

void nlSolverParameteri(NLenum pname, NLint param) {
    nlCheckState(NL_STATE_INITIAL) ;
    switch(pname) {
    case NL_SOLVER: {
        nlCurrentContext->solver = (NLenum)param ;
    } break ;
    case NL_NB_VARIABLES: {
        nl_assert(param > 0) ;
        nlCurrentContext->nb_variables = (NLuint)param ;
    } break ;
    case NL_LEAST_SQUARES: {
        nlCurrentContext->least_squares = (NLboolean)param ;
    } break ;
    case NL_MAX_ITERATIONS: {
        nl_assert(param > 0) ;
        nlCurrentContext->max_iterations = (NLuint)param ;
    } break ;
    case NL_THRESHOLD: {
        nl_assert(param >= 0) ;
        nlCurrentContext->threshold = (NLdouble)param ;
    } break ;
    case NL_OMEGA: {
        nl_range_assert(param,1,2) ;
        nlCurrentContext->omega = (NLdouble)param ;
    } break ;
    case NL_SYMMETRIC: {
        nlCurrentContext->symmetric = (NLboolean)param ;        
    }
    case NL_INNER_ITERATIONS: {
        nl_assert(param > 0) ;
        nlCurrentContext->inner_iterations = (NLuint)param ;
    } break ;
    case NL_PRECONDITIONER: {
        nlCurrentContext->preconditioner = (NLuint)param ;        
    } break ;
    default: {
        nl_assert_not_reached ;
    } break ;
    }
}

void nlRowParameterd(NLenum pname, NLdouble param) {
    nlCheckState(NL_STATE_MATRIX) ;
    switch(pname) {
    case NL_RIGHT_HAND_SIDE: {
        nlCurrentContext->right_hand_side = param ;
    } break ;
    case NL_ROW_SCALING: {
        nlCurrentContext->row_scaling = param ;
    } break ;
    }
}

void nlRowParameteri(NLenum pname, NLint param) {
    nlCheckState(NL_STATE_MATRIX) ;
    switch(pname) {
    case NL_RIGHT_HAND_SIDE: {
        nlCurrentContext->right_hand_side = (NLdouble)param ;
    } break ;
    case NL_ROW_SCALING: {
        nlCurrentContext->row_scaling = (NLdouble)param ;
    } break ;
    }
}

void nlGetBooleanv(NLenum pname, NLboolean* params) {
    switch(pname) {
    case NL_LEAST_SQUARES: {
        *params = nlCurrentContext->least_squares ;
    } break ;
    case NL_SYMMETRIC: {
        *params = nlCurrentContext->symmetric ;
    } break ;
    default: {
        nl_assert_not_reached ;
    } break ;
    }
}

void nlGetDoublev(NLenum pname, NLdouble* params) {
    switch(pname) {
    case NL_SOLVER: {
        *params = (NLdouble)(nlCurrentContext->solver) ;
    } break ;
    case NL_NB_VARIABLES: {
        *params = (NLdouble)(nlCurrentContext->nb_variables) ;
    } break ;
    case NL_LEAST_SQUARES: {
        *params = (NLdouble)(nlCurrentContext->least_squares) ;
    } break ;
    case NL_MAX_ITERATIONS: {
        *params = (NLdouble)(nlCurrentContext->max_iterations) ;
    } break ;
    case NL_THRESHOLD: {
        *params = (NLdouble)(nlCurrentContext->threshold) ;
    } break ;
    case NL_OMEGA: {
        *params = (NLdouble)(nlCurrentContext->omega) ;
    } break ;
    case NL_SYMMETRIC: {
        *params = (NLdouble)(nlCurrentContext->symmetric) ;
    } break ;
    case NL_USED_ITERATIONS: {
        *params = (NLdouble)(nlCurrentContext->used_iterations) ;
    } break ;
    case NL_ERROR: {
        *params = (NLdouble)(nlCurrentContext->error) ;
    } break ;
    case NL_ELAPSED_TIME: {
        *params = (NLdouble)(nlCurrentContext->elapsed_time) ;        
    } break ;
    case NL_PRECONDITIONER: {
        *params = (NLdouble)(nlCurrentContext->preconditioner) ;        
    } break ;
    default: {
        nl_assert_not_reached ;
    } break ;
    }
}

void nlGetIntergerv(NLenum pname, NLint* params) {
    switch(pname) {
    case NL_SOLVER: {
        *params = (NLint)(nlCurrentContext->solver) ;
    } break ;
    case NL_NB_VARIABLES: {
        *params = (NLint)(nlCurrentContext->nb_variables) ;
    } break ;
    case NL_LEAST_SQUARES: {
        *params = (NLint)(nlCurrentContext->least_squares) ;
    } break ;
    case NL_MAX_ITERATIONS: {
        *params = (NLint)(nlCurrentContext->max_iterations) ;
    } break ;
    case NL_THRESHOLD: {
        *params = (NLint)(nlCurrentContext->threshold) ;
    } break ;
    case NL_OMEGA: {
        *params = (NLint)(nlCurrentContext->omega) ;
    } break ;
    case NL_SYMMETRIC: {
        *params = (NLint)(nlCurrentContext->symmetric) ;
    } break ;
    case NL_USED_ITERATIONS: {
        *params = (NLint)(nlCurrentContext->used_iterations) ;
    } break ;
    case NL_PRECONDITIONER: {
        *params = (NLint)(nlCurrentContext->preconditioner) ;        
    } break ;
    default: {
        nl_assert_not_reached ;
    } break ;
    }
}

/************************************************************************************/
/* Enable / Disable */

void nlEnable(NLenum pname) {
    switch(pname) {
    case NL_NORMALIZE_ROWS: {
        nl_assert(nlCurrentContext->state != NL_STATE_ROW) ;
        nlCurrentContext->normalize_rows = NL_TRUE ;
    } break ;
    default: {
        nl_assert_not_reached ;
    }
    }
}

void nlDisable(NLenum pname) {
    switch(pname) {
    case NL_NORMALIZE_ROWS: {
        nl_assert(nlCurrentContext->state != NL_STATE_ROW) ;
        nlCurrentContext->normalize_rows = NL_FALSE ;
    } break ;
    default: {
        nl_assert_not_reached ;
    }
    }
}

NLboolean nlIsEnabled(NLenum pname) {
    switch(pname) {
    case NL_NORMALIZE_ROWS: {
        return nlCurrentContext->normalize_rows ;
    } break ;
    default: {
        nl_assert_not_reached ;
    }
    }
    return NL_FALSE ;
}

/************************************************************************************/
/* NL functions */

void  nlSetFunction(NLenum pname, void* param) {
    switch(pname) {
    case NL_FUNC_SOLVER:
        nlCurrentContext->solver_func = (NLSolverFunc)(param)  ;
        break ;
    case NL_FUNC_MATRIX:
        nlCurrentContext->matrix_vector_prod = (NLMatrixFunc)(param)  ;
        nlCurrentContext->solver = NL_SOLVER_USER ;
        break ;
    case NL_FUNC_PRECONDITIONER:
        nlCurrentContext->precond_vector_prod = (NLMatrixFunc)(param)  ;
        nlCurrentContext->preconditioner = NL_PRECOND_USER ;
        break ;
    default:
        nl_assert_not_reached ;
    }
}

void nlGetFunction(NLenum pname, void** param) {
    switch(pname) {
    case NL_FUNC_SOLVER:
        *param = nlCurrentContext->solver_func ;
        break ;
    case NL_FUNC_MATRIX:
        *param = nlCurrentContext->matrix_vector_prod ;
        break ;
    case NL_FUNC_PRECONDITIONER:
        *param = nlCurrentContext->precond_vector_prod ;
        break ;
    default:
        nl_assert_not_reached ;
    }
}

/************************************************************************************/
/* Get/Set Lock/Unlock variables */

void nlSetVariable(NLuint index, NLdouble value) {
    nlCheckState(NL_STATE_SYSTEM) ;
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1) ;
    nlCurrentContext->variable[index].value = value ;    
}

NLdouble nlGetVariable(NLuint index) {
    nl_assert(nlCurrentContext->state != NL_STATE_INITIAL) ;
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1) ;
    return nlCurrentContext->variable[index].value ;
}

void nlLockVariable(NLuint index) {
    nlCheckState(NL_STATE_SYSTEM) ;
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1) ;
    nlCurrentContext->variable[index].locked = NL_TRUE ;
}

void nlUnlockVariable(NLuint index) {
    nlCheckState(NL_STATE_SYSTEM) ;
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1) ;
    nlCurrentContext->variable[index].locked = NL_FALSE ;
}

NLboolean nlVariableIsLocked(NLuint index) {
    nl_assert(nlCurrentContext->state != NL_STATE_INITIAL) ;
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1) ;
    return nlCurrentContext->variable[index].locked  ;
}

/************************************************************************************/
/* System construction */

void nlVariablesToVector() {
    NLuint i ;
    nl_assert(nlCurrentContext->alloc_x) ;
    nl_assert(nlCurrentContext->alloc_variable) ;
    for(i=0; i<nlCurrentContext->nb_variables; i++) {
        NLVariable* v = &(nlCurrentContext->variable[i]) ;
        if(!v->locked) {
            nl_assert(v->index < nlCurrentContext->n) ;
            nlCurrentContext->x[v->index] = v->value ;
        }
    }
}

void nlVectorToVariables() {
    NLuint i ;
    nl_assert(nlCurrentContext->alloc_x) ;
    nl_assert(nlCurrentContext->alloc_variable) ;
    for(i=0; i<nlCurrentContext->nb_variables; i++) {
        NLVariable* v = &(nlCurrentContext->variable[i]) ;
        if(!v->locked) {
            nl_assert(v->index < nlCurrentContext->n) ;
            v->value = nlCurrentContext->x[v->index] ;
        }
    }
}


void nlBeginSystem() {
    nlTransition(NL_STATE_INITIAL, NL_STATE_SYSTEM) ;
    nl_assert(nlCurrentContext->nb_variables > 0) ;
    nlCurrentContext->variable = NL_NEW_ARRAY(
        NLVariable, nlCurrentContext->nb_variables
    ) ;
    nlCurrentContext->alloc_variable = NL_TRUE ;
}

void nlEndSystem() {
    nlTransition(NL_STATE_MATRIX_CONSTRUCTED, NL_STATE_SYSTEM_CONSTRUCTED) ;    
}

void nlBeginMatrix() {
    NLuint i ;
    NLuint n = 0 ;
    NLenum storage = NL_MATRIX_STORE_ROWS ;

    nlTransition(NL_STATE_SYSTEM, NL_STATE_MATRIX) ;

    for(i=0; i<nlCurrentContext->nb_variables; i++) {
        if(!nlCurrentContext->variable[i].locked) {
            nlCurrentContext->variable[i].index = n ;
            n++ ;
        } else {
            nlCurrentContext->variable[i].index = ~0 ;
        }
    }

    nlCurrentContext->n = n ;

    /* SSOR preconditioner requires rows and columns */
    if(nlCurrentContext->preconditioner == NL_PRECOND_SSOR) {
        storage = (storage | NL_MATRIX_STORE_COLUMNS) ;
    }

    /* a least squares problem results in a symmetric matrix */
    if(nlCurrentContext->least_squares 
        && !nlSolverIsCNC(nlCurrentContext->solver)) {
        nlCurrentContext->symmetric = NL_TRUE ;
    }

    if(nlCurrentContext->symmetric) {
        storage = (storage | NL_MATRIX_STORE_SYMMETRIC) ;
    }

    /* SuperLU storage does not support symmetric storage */
    if(
        nlCurrentContext->solver == NL_SUPERLU_EXT       ||
        nlCurrentContext->solver == NL_PERM_SUPERLU_EXT  ||
        nlCurrentContext->solver == NL_SYMMETRIC_SUPERLU_EXT 
    ) {
        storage = (storage & ~NL_SYMMETRIC) ;
    }

    nlSparseMatrixConstruct(&nlCurrentContext->M, n, n, storage) ;
    nlCurrentContext->alloc_M = NL_TRUE ;

    nlCurrentContext->x = NL_NEW_ARRAY(NLdouble, n) ;
    nlCurrentContext->alloc_x = NL_TRUE ;
    
    nlCurrentContext->b = NL_NEW_ARRAY(NLdouble, n) ;
    nlCurrentContext->alloc_b = NL_TRUE ;

    nlVariablesToVector() ;

    nlRowColumnConstruct(&nlCurrentContext->af) ;
    nlCurrentContext->alloc_af = NL_TRUE ;
    nlRowColumnConstruct(&nlCurrentContext->al) ;
    nlCurrentContext->alloc_al = NL_TRUE ;
    nlRowColumnConstruct(&nlCurrentContext->xl) ;
    nlCurrentContext->alloc_xl = NL_TRUE ;

    nlCurrentContext->current_row = 0 ;
}

void nlEndMatrix() {
    nlTransition(NL_STATE_MATRIX, NL_STATE_MATRIX_CONSTRUCTED) ;    
    
    nlRowColumnDestroy(&nlCurrentContext->af) ;
    nlCurrentContext->alloc_af = NL_FALSE ;
    nlRowColumnDestroy(&nlCurrentContext->al) ;
    nlCurrentContext->alloc_al = NL_FALSE ;
    nlRowColumnDestroy(&nlCurrentContext->xl) ;
    nlCurrentContext->alloc_al = NL_FALSE ;
    
    if(!nlCurrentContext->least_squares) {
        nl_assert(
            nlCurrentContext->current_row == 
            nlCurrentContext->n
        ) ;
    }
}

void nlBeginRow() {
    nlTransition(NL_STATE_MATRIX, NL_STATE_ROW) ;
    nlRowColumnZero(&nlCurrentContext->af) ;
    nlRowColumnZero(&nlCurrentContext->al) ;
    nlRowColumnZero(&nlCurrentContext->xl) ;
}

void nlScaleRow(NLdouble s) {
    NLRowColumn*    af = &nlCurrentContext->af ;
    NLRowColumn*    al = &nlCurrentContext->al ;
    NLuint nf            = af->size ;
    NLuint nl            = al->size ;
    NLuint i ;
    for(i=0; i<nf; i++) {
        af->coeff[i].value *= s ;
    }
    for(i=0; i<nl; i++) {
        al->coeff[i].value *= s ;
    }
    nlCurrentContext->right_hand_side *= s ;
}

void nlNormalizeRow(NLdouble weight) {
    NLRowColumn*    af = &nlCurrentContext->af ;
    NLRowColumn*    al = &nlCurrentContext->al ;
    NLuint nf            = af->size ;
    NLuint nl            = al->size ;
    NLuint i ;
    NLdouble norm = 0.0 ;
    for(i=0; i<nf; i++) {
        norm += af->coeff[i].value * af->coeff[i].value ;
    }
    for(i=0; i<nl; i++) {
        norm += al->coeff[i].value * al->coeff[i].value ;
    }
    norm = sqrt(norm) ;
    nlScaleRow(weight / norm) ;
}

void nlEndRow() {
    NLRowColumn*    af = &nlCurrentContext->af ;
    NLRowColumn*    al = &nlCurrentContext->al ;
    NLRowColumn*    xl = &nlCurrentContext->xl ;
    NLSparseMatrix* M  = &nlCurrentContext->M  ;
    NLdouble* b        = nlCurrentContext->b ;
    NLuint nf          = af->size ;
    NLuint nl          = al->size ;
    NLuint current_row = nlCurrentContext->current_row ;
    NLuint i ;
    NLuint j ;
    NLdouble S ;
    nlTransition(NL_STATE_ROW, NL_STATE_MATRIX) ;

    if(nlCurrentContext->normalize_rows) {
        nlNormalizeRow(nlCurrentContext->row_scaling) ;
    } else {
        nlScaleRow(nlCurrentContext->row_scaling) ;
    }
    // if least_squares : we want to solve
    // A'A x = A'b
    //
    if(nlCurrentContext->least_squares) {
        for(i=0; i<nf; i++) {
            for(j=0; j<nf; j++) {
                nlSparseMatrixAdd(
                    M, af->coeff[i].index, af->coeff[j].index,
                    af->coeff[i].value * af->coeff[j].value
                ) ;
            }
        }
        S = -nlCurrentContext->right_hand_side ;
        for(j=0; j<nl; j++) {
            S += al->coeff[j].value * xl->coeff[j].value ;
        }
        for(i=0; i<nf; i++) {
            b[ af->coeff[i].index ] -= af->coeff[i].value * S ;
        }
    } else {
        for(i=0; i<nf; i++) {
            nlSparseMatrixAdd(
                M, current_row, af->coeff[i].index, af->coeff[i].value
            ) ;
        }
        b[current_row] = -nlCurrentContext->right_hand_side ;
        for(i=0; i<nl; i++) {
            b[current_row] -= al->coeff[i].value * xl->coeff[i].value ;
        }
    }
    nlCurrentContext->current_row++ ;
    nlCurrentContext->right_hand_side = 0.0 ;    
    nlCurrentContext->row_scaling     = 1.0 ;
}

void nlCoefficient(NLuint index, NLdouble value) {
    NLVariable* v = NULL ;
    nlCheckState(NL_STATE_ROW) ;
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1) ;
    v = &(nlCurrentContext->variable[index]) ;
    if(v->locked) {
        nlRowColumnAppend(&(nlCurrentContext->al), 0, value) ;
        nlRowColumnAppend(&(nlCurrentContext->xl), 0, v->value) ;
    } else {
        nlRowColumnAppend(&(nlCurrentContext->af), v->index, value) ;
    }
}

void nlBegin(NLenum prim) {
    switch(prim) {
    case NL_SYSTEM: {
        nlBeginSystem() ;
    } break ;
    case NL_MATRIX: {
        nlBeginMatrix() ;
    } break ;
    case NL_ROW: {
        nlBeginRow() ;
    } break ;
    default: {
        nl_assert_not_reached ;
    }
    }
}

void nlEnd(NLenum prim) {
    switch(prim) {
    case NL_SYSTEM: {
        nlEndSystem() ;
    } break ;
    case NL_MATRIX: {
        nlEndMatrix() ;
    } break ;
    case NL_ROW: {
        nlEndRow() ;
    } break ;
    default: {
        nl_assert_not_reached ;
    }
    }
}

/************************************************************************/
/* nlSolve() driver routine */

NLboolean nlSolve() {
    NLboolean result ;
    NLdouble start_time = nlCurrentTime() ; 
    nlCheckState(NL_STATE_SYSTEM_CONSTRUCTED) ;
    nlCurrentContext->elapsed_time = 0 ;
    result =  nlCurrentContext->solver_func() ;
    nlVectorToVariables() ;
    nlCurrentContext->elapsed_time = nlCurrentTime() - start_time ;
    nlTransition(NL_STATE_SYSTEM_CONSTRUCTED, NL_STATE_SOLVED) ;
    return result ;
}

