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

#include <NL/nl_cnc_gpu_cuda.h>
#include <NL/nl_context.h>

NLboolean nlSolverIsCNC(NLint solver){
    return solver == NL_CNC_FLOAT_CRS 
        || solver == NL_CNC_DOUBLE_CRS
        || solver == NL_CNC_FLOAT_BCRS2
        || solver == NL_CNC_DOUBLE_BCRS2        
        || solver == NL_CNC_FLOAT_ELL        
        || solver == NL_CNC_DOUBLE_ELL        
        || solver == NL_CNC_FLOAT_HYB        
        || solver == NL_CNC_DOUBLE_HYB;        
}


/************************************************************************/
/* CNC wrapper */

#ifdef NL_USE_CNC

NLuint nlSolve_CNC() {
    unsigned int i;
    NLdouble* b        = nlCurrentContext->b ;
    NLdouble* x        = nlCurrentContext->x ;
    NLdouble  eps      = nlCurrentContext->threshold ;
    NLuint    max_iter = nlCurrentContext->max_iterations ;
    NLSparseMatrix *M  = &(nlCurrentContext->M);
    
    // local variables for the final error computation
    NLuint val_ret;
    NLdouble * Ax=NL_NEW_ARRAY(NLdouble,nlCurrentContext->n);
    NLdouble accu     = 0.0;
    NLdouble b_square = 0.0;
    //nl_assert( M->n == nlCurrentContext->n);  
    
    // call to cnc solver
    val_ret=cnc_solve_cg(M, b, x, max_iter, eps, nlCurrentContext->solver);
    
    // compute the final error 
    nlCurrentContext->matrix_vector_prod(x,Ax);
    for(i = 0 ; i < M->n ; ++i) { 
        accu     +=(Ax[i]-b[i])*(Ax[i]-b[i]);
        b_square += b[i]*b[i]; 
    } 
    printf("in OpenNL : ||Ax-b||/||b|| = %e\n",sqrt(accu)/sqrt(b_square));
    // cleaning 
    NL_DELETE_ARRAY(Ax);
    return val_ret;
}

#else

NLuint nlSolve_CNC() {
    nl_assert_not_reached ;
}

#endif
