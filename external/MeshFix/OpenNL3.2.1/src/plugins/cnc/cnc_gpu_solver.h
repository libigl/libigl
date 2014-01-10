/*
 *  Copyright (c) 2004-2010, Bruno Levy
 *  All rights reserved.
 *
 *  CNC: Concurrent Number Cruncher, original code by Luc Buatois
 *  Copyright (C) 2008-2010 GOCAD/ASGA, INRIA/ALICE
 *
 *  Sparse matrix-vector multiplication (SpMV) CUDA kernels based on code
 *  by Nathan Bell and Michael Garland at NVIDIA.
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

#ifndef CNC_GPU_SOLVER_H
#define CNC_GPU_SOLVER_H

#include <NL/nl.h>
#include <NL/nl_cnc_gpu_cuda.h>
#ifdef __cplusplus
#include "cnc_utils.h"
#include "cnc_cublas_utils.h"
#include "cnc_sparse_matrix_ell.h"
#include "cnc_sparse_matrix_coo.h"
#include "cnc_sparse_matrix_crs.h"
#include "cnc_sparse_matrix_bcrs.h"
#include "cnc_sparse_matrix_hyb.h"
#include "cnc_timer.h"

#include <vector>
#include <algorithm>



//------------------------------------------------------------------------------//
// Implementation of the Conjugate Gradient Solver on the GPU                   //
// Matrix \in { CNCSparseMatrixCRS<float>, CNCSparseMatrixCRS<double>,          //
//             CNCSparseMatrixBCRS<float,2,2>, CNCSparseMatrixBCRS<double,2,2>  //
//             CNCSparseMatrixELL<float>, CNCSparseMatrixELL<double>            //
//             CNCSparseMatrixCOO<float>, CNCSparseMatrixCOO<double>            //
//             CNCSparseMatrixHYB<float>, CNCSparseMatrixHYB<double>            //
// TypeVector \in { float, double }                                             //
//------------------------------------------------------------------------------//
template<typename Matrix,typename TypeVector> 
inline NLuint solve_cg_internal (
		Matrix &A,
		CNCArray1d<TypeVector>  & b,
		CNCArray1d<TypeVector>  & x,
		const unsigned int nb_iter_max,
		const double epsilon) {


	CNCTimer time ;
	double full_time = time.GetTime() ;

	unsigned int N = x.size();
	// init and check cublas
	cublasStatus st = cublasInit () ;
	CNCgetError ( st ) ;


	// vars to be defined specifically for each storage format (CRS/BCRS...)
	const typename Matrix::type_value * diag_matrix = A.diagonal() ;

	TypeVector *gpu_r = NULL ;
	TypeVector *gpu_d = NULL ;
	TypeVector *gpu_h = NULL ;
	TypeVector *gpu_Ad= NULL ;
	TypeVector *gpu_diag_inv = NULL ;
	TypeVector *gpu_b = NULL ;
	TypeVector *gpu_x = NULL ;

	// matrix allocation and upload
	A.gpu_allocate_and_upload() ;

	new_device_vector ( &gpu_r, N+16 );
	new_device_vector ( &gpu_d, N+16 );
	new_device_vector ( &gpu_h, N+16 );
	new_device_vector ( &gpu_Ad, N+16 );
	new_device_vector ( &gpu_diag_inv, N+16 );
	new_device_vector ( &gpu_x, N+16 );
	new_device_vector ( &gpu_b, N+16 );

	// building the Jacobi preconditionner
 
	CNCArray1d<TypeVector>  cpu_diag_inv ( N+16 ) ;
	for(unsigned int i=0; i<N; i++) {
		cpu_diag_inv[i] = static_cast<TypeVector>(
		    (i >= N || diag_matrix[i] == 0.0) ? 1.0 : 1.0 / diag_matrix[i]
		    ) ;
	}

	copy_host_to_device_vector ( cpu_diag_inv, gpu_diag_inv );
	copy_host_to_device_vector ( x           , gpu_x        );
	copy_host_to_device_vector ( b           , gpu_b        );

	unsigned int its=0;
	TypeVector alpha, beta;

	printf ( "------------------------------------------------------------\n" ) ;
	printf ( "End init of CG solver\n" ) ;
	printf ( "------------------------------------------------------------\n" ) ;
	printf ( "Start of CG solver\n" ) ;
	printf ( "------------------------------------------------------------\n" ) ;
	
	// r = A*x
	A.gpu_mult ( gpu_x, gpu_r, N ) ;

	// r = b - A*x
	wrapXaxpy ( N, -1.0, gpu_b, 1, gpu_r, 1 ) ;
	wrapXscal ( N, -1.0, gpu_r, 1 ) ;

        
	const unsigned int BLOCK_SIZE = 256;
    const dim3 grid = make_large_grid(N, BLOCK_SIZE);
	// d = M^(-1) * r
    // Hadamard product between diagonal preconditioner and r
	CNCVecVecMultKernel<<<grid, BLOCK_SIZE>>>( N, gpu_diag_inv, gpu_r, gpu_d );
	// cur_err = rTypeVector*d err = epsilon^2 * rTypeVector*d
	TypeVector cur_err = wrapXdot ( N, gpu_r, 1, gpu_d, 1 ) ;

    // using same stopping criterion than NL_CG
    #define  TRAD_ERROR  
        
	#ifndef  TRAD_ERROR  
    TypeVector err = static_cast<TypeVector>(cur_err * epsilon * epsilon) ;
	while ( cur_err > err && its < nb_iter_max) {
		if(!(its & 31)) {
			printf ( "%d : %.10e -- %.10e\n", its, cur_err, err ) ;
		}
        #else
	// cur_err = ||r||^2 and err= epsilon^2 ||b|| 
	TypeVector trad_err = wrapXdot ( N, gpu_r, 1, gpu_r, 1 ) ;
	
	TypeVector b_square = wrapXdot ( N, gpu_b, 1, gpu_b, 1 ) ;
        
	TypeVector err = static_cast<TypeVector>(b_square * epsilon * epsilon) ;
	
	while ( trad_err > err && its < nb_iter_max) {
		if(!(its & 31)) {
			printf ( "%d : %.10e -- %.10e\n", its, trad_err, err ) ;
		}
        #endif        

		// Ad = A*d
		A.gpu_mult ( gpu_d, gpu_Ad, N) ;

		// alpha = cur_err / (dTypeVector*Ad)
		alpha = cur_err / wrapXdot ( N, gpu_d, 1, gpu_Ad, 1 ) ;

		// x = x + alpha * d
		wrapXaxpy ( N, alpha, gpu_d, 1, gpu_x, 1 ) ;

		// r = r - alpha * Ad
		wrapXaxpy ( N, -alpha, gpu_Ad, 1, gpu_r, 1 ) ;

		// h = M^(-1) * r
        // Hadamard product between diagonal preconditioner and r
        CNCVecVecMultKernel<<<grid, BLOCK_SIZE>>>( N, gpu_diag_inv, gpu_r, gpu_h );

		TypeVector old_err = cur_err ;

		// cur_err = rT * h
		cur_err = wrapXdot ( N, gpu_r, 1, gpu_h, 1 ) ;

		beta = cur_err / old_err ;

		// d = h + beta * d
		wrapXscal ( N, beta, gpu_d, 1 ) ;
		
		wrapXaxpy ( N, 1.0, gpu_h, 1, gpu_d, 1 ) ;
		++its;

        // using same stopping criterion than NL_CG
        #ifdef  TRAD_ERROR
        trad_err = wrapXdot (N, gpu_r, 1, gpu_r, 1 ) ;
        #endif

	}

	printf ( "------------------------------------------------------------\n" ) ;
	if ( its==nb_iter_max ) {
		printf ( "Maximum #itr reached: SOLVER DID NOT CONVERGE !!!\n" ) ; 
		printf ( "------------------------------------------------------------\n" ) ;
	}
	
	double seconds = time.GetElapsedTime(full_time) ;
    double sec_per_iteration = seconds / its;
    double GFLOPs = (sec_per_iteration == 0) ? 0 : 
        ((2.0 * (double) (A.data_size())+(double)(11*N))/sec_per_iteration)/1e9;
    printf("Nathan Bell GFLOPs %2.2f GFLOP/s  Time : %g\n",GFLOPs, seconds);

	copy_device_to_host_vector( gpu_x,  x ) ;

    // check test : value of ||Ax-b|| on CPU
	CNCArray1d<TypeVector> Ax ( N );
    A.mult( x, Ax);
    TypeVector accu = static_cast<TypeVector>(0);
	for(unsigned int i = 0; i < N ; ++i){
        accu += (Ax[i]-b[i])*(Ax[i]-b[i]);
    }
	printf("in CNC    : ||Ax-b||/||b|| = %e\n",sqrt_gen(accu)/sqrt_gen(b_square));
    
    
    delete_device_vector ( &gpu_r ) ;
	delete_device_vector ( &gpu_d ) ;
	delete_device_vector ( &gpu_h  ) ;
	delete_device_vector ( &gpu_Ad ) ;
	delete_device_vector ( &gpu_diag_inv ) ;
	delete_device_vector ( &gpu_x ) ;
	delete_device_vector ( &gpu_b ) ;

	A.gpu_deallocate () ;

	st = cublasShutdown() ;
	CNCgetError ( st ) ;

	return its ;
}

extern "C" {
#endif
NLuint cnc_solve_cg (NLSparseMatrix *A, NLdouble *b, NLdouble *x,
        NLuint nb_iter_max, NLdouble epsilon, NLint solver_type );

#ifdef __cplusplus
}
#endif

#endif

//---------------------------------------------------------------------------//
