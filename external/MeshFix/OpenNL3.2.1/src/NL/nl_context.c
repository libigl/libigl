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

#include <NL/nl_context.h>
#include <NL/nl_iterative_solvers.h>
#include <NL/nl_preconditioners.h>
#include <NL/nl_superlu.h>
#include <NL/nl_cnc_gpu_cuda.h>

NLContextStruct* nlCurrentContext = NULL ;

void nlMatrixVectorProd_default(NLdouble* x, NLdouble* y) {
    nlSparseMatrixMult(&(nlCurrentContext->M), x, y) ;
}

NLContext nlNewContext() {
    NLContextStruct* result    = NL_NEW(NLContextStruct) ;
    result->state              = NL_STATE_INITIAL ;
    result->solver             = NL_BICGSTAB ;
    result->max_iterations     = 100 ;
    result->threshold          = 1e-6 ;
    result->omega              = 1.5 ;
    result->row_scaling        = 1.0 ;
    result->right_hand_side    = 0.0 ;
    result->inner_iterations   = 5 ;
    result->matrix_vector_prod = nlMatrixVectorProd_default ;
    result->solver_func        = nlDefaultSolver ;
    nlMakeCurrent(result) ;
    return result ;
}

void nlDeleteContext(NLContext context_in) {
    NLContextStruct* context = (NLContextStruct*)(context_in) ;
    if(nlCurrentContext == context) {
        nlCurrentContext = NULL ;
    }
    if(context->alloc_M) {
        nlSparseMatrixDestroy(&context->M) ;
    }
    if(context->alloc_af) {
        nlRowColumnDestroy(&context->af) ;
    }
    if(context->alloc_al) {
        nlRowColumnDestroy(&context->al) ;
    }
    if(context->alloc_xl) {
        nlRowColumnDestroy(&context->xl) ;
    }
    if(context->alloc_variable) {
        NL_DELETE_ARRAY(context->variable) ;
    }
    if(context->alloc_x) {
        NL_DELETE_ARRAY(context->x) ;
    }
    if(context->alloc_b) {
        NL_DELETE_ARRAY(context->b) ;
    }

#ifdef NL_PARANOID
    NL_CLEAR(NLContextStruct, context) ;
#endif
    NL_DELETE(context) ;
}

void nlMakeCurrent(NLContext context) {
    nlCurrentContext = (NLContextStruct*)(context) ;
}

NLContext nlGetCurrent() {
    return nlCurrentContext ;
}

/************************************************************************/
/* Finite state automaton   */

void nlCheckState(NLenum state) {
    nl_assert(nlCurrentContext->state == state) ;
}

void nlTransition(NLenum from_state, NLenum to_state) {
    nlCheckState(from_state) ;
    nlCurrentContext->state = to_state ;
}

/************************************************************************/
/* Default solver */

static void nlSetupPreconditioner() {
    switch(nlCurrentContext->preconditioner) {
    case NL_PRECOND_NONE:
        nlCurrentContext->precond_vector_prod = NULL ;
        break ;
    case NL_PRECOND_JACOBI:
        nlCurrentContext->precond_vector_prod = nlPreconditioner_Jacobi ;
        break ;
    case NL_PRECOND_SSOR:
        nlCurrentContext->precond_vector_prod = nlPreconditioner_SSOR ;
        break ;
    default:
        nl_assert_not_reached ;
        break ;
    }
    /* Check compatibility between solver and preconditioner */
    if(
        nlCurrentContext->solver == NL_BICGSTAB && 
        nlCurrentContext->preconditioner == NL_PRECOND_SSOR
    ) {
        nlWarning(
            "nlSolve", 
            "cannot use SSOR preconditioner with non-symmetric matrix, switching to Jacobi"
        ) ;
        nlCurrentContext->preconditioner = NL_PRECOND_JACOBI ;        
        nlCurrentContext->precond_vector_prod = nlPreconditioner_Jacobi ;
    }
    if(
        nlCurrentContext->solver == NL_GMRES && 
        nlCurrentContext->preconditioner != NL_PRECOND_NONE
    ) {
        nlWarning("nlSolve", "Preconditioner not implemented yet for GMRES") ;
        nlCurrentContext->preconditioner = NL_PRECOND_NONE ;        
        nlCurrentContext->precond_vector_prod = NULL ;
    }
    if(
        nlCurrentContext->solver == NL_SUPERLU_EXT && 
        nlCurrentContext->preconditioner != NL_PRECOND_NONE
    ) {
        nlWarning("nlSolve", "Preconditioner not implemented yet for SUPERLU") ;
        nlCurrentContext->preconditioner = NL_PRECOND_NONE ;        
        nlCurrentContext->precond_vector_prod = NULL ;
    }
    if(
        nlCurrentContext->solver == NL_PERM_SUPERLU_EXT && 
        nlCurrentContext->preconditioner != NL_PRECOND_NONE
    ) {
        nlWarning("nlSolve", "Preconditioner not implemented yet for PERMSUPERLU") ;
        nlCurrentContext->preconditioner = NL_PRECOND_NONE ;        
        nlCurrentContext->precond_vector_prod = NULL ;
    }
    if(
        nlCurrentContext->solver == NL_SYMMETRIC_SUPERLU_EXT && 
        nlCurrentContext->preconditioner != NL_PRECOND_NONE
    ) {
        nlWarning("nlSolve", "Preconditioner not implemented yet for PERMSUPERLU") ;
        nlCurrentContext->preconditioner = NL_PRECOND_NONE ;        
        nlCurrentContext->precond_vector_prod = NULL ;
    }
}

NLboolean nlDefaultSolver() {
    NLboolean result = NL_TRUE ;
    nlSetupPreconditioner() ;
    switch(nlCurrentContext->solver) {
    case NL_CG: {
        if(nlCurrentContext->preconditioner == NL_PRECOND_NONE) {
            nlCurrentContext->used_iterations = nlSolve_CG() ;
        } else {
            nlCurrentContext->used_iterations = nlSolve_CG_precond() ;
        }
    } break ;
    case NL_BICGSTAB: {
        if(nlCurrentContext->preconditioner == NL_PRECOND_NONE) {
            nlCurrentContext->used_iterations = nlSolve_BICGSTAB() ;
        } else {
            nlCurrentContext->used_iterations = nlSolve_BICGSTAB_precond() ;
        }
    } break ;
    case NL_GMRES: {
        nlCurrentContext->used_iterations = nlSolve_GMRES() ;
    } break ;
    case NL_CNC_FLOAT_CRS:
    case NL_CNC_DOUBLE_CRS:
    case NL_CNC_FLOAT_BCRS2:
    case NL_CNC_DOUBLE_BCRS2:
    case NL_CNC_FLOAT_ELL:
    case NL_CNC_DOUBLE_ELL:
    case NL_CNC_FLOAT_HYB:
    case NL_CNC_DOUBLE_HYB: {
        nlCurrentContext->used_iterations = nlSolve_CNC() ;
    } break ;
    case NL_SUPERLU_EXT: 
    case NL_PERM_SUPERLU_EXT: 
    case NL_SYMMETRIC_SUPERLU_EXT: {
        result = nlSolve_SUPERLU() ;
    } break ;
    default:
        nl_assert_not_reached ;
    }
    return result ;
}
