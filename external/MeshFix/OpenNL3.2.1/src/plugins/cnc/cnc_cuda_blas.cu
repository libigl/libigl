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

// Header files
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cnc_kernels.h"
#include "cnc_gpu_solver.h"
#include "cnc_texture.h"




// a convert fonctor to use with the transform algorithm
// when copying from double vector to float vector and vice-versa
template <typename Arg,typename Result>
struct convert : std::unary_function<Arg,Result> {
    Result operator() (const Arg x) {
        return  static_cast<Result>(x);
    }   
}; 


//---------------------------------------------------------------------------//
// kernel for matrix-product vector with a BCRS<float,2,2> matrix            //  
//---------------------------------------------------------------------------//
__global__ void CNCMat2x2VecMult4Kernel (
        float4 * matrix,
        unsigned int size_matrix,
        uint2 * rowptr,
        unsigned int size_rowptr,
        unsigned int * colind,
        unsigned int size_colind,
        const float2 * x,
        float2 * b,
        unsigned int size_vec ) {

	// Thread index
	const unsigned int index = large_grid_thread_id(void) ;

	if ( index<<1 < size_vec ) {
		uint2 rowptr_bounds = rowptr[index] ;
		float2 res ;
		res.x = res.y = 0.0f ;

		unsigned int ci = 0 ;
		float2 x_vec ;
		float4 mat_vec ;

		// for each block of the block_row, mult
		for ( int i=rowptr_bounds.x; i<rowptr_bounds.y; i++ ) { 
			ci = colind[i] ;
			mat_vec = matrix[i] ;
			#ifndef USE_TEXTURE
			x_vec = x[ci] ;
			#else
			x_vec = fetch_x(ci,x) ;
			#endif

			res.x += mat_vec.x*x_vec.x+mat_vec.y*x_vec.y ;
			res.y += mat_vec.z*x_vec.x+mat_vec.w*x_vec.y ;
		}
		b[index] = res ;
	}
}


//---------------------------------------------------------------------------//
// kernel for matrix-product vector with a BCRS<double,2,2> matrix           //  
//---------------------------------------------------------------------------//

__global__ void CNCMat2x2VecMult4Kernel (
        double2 * mat0,
        double2 *mat1 ,
        unsigned int size_matrix,
        uint2 * rowptr,
        unsigned int size_rowptr,
        unsigned int * colind,
        unsigned int size_colind,
        const double2 * x,
        double2 * b,
        unsigned int size_vec ) {

	// Thread index
	const unsigned int index = large_grid_thread_id(void); 

	if ( index<<1 < size_vec ) {

		uint2 rowptr_bounds = rowptr[index] ;
		double2 res ;
		res.x = res.y = 0.0;

		unsigned int ci = 0 ;
		double2 row0;
		double2 row1;
		double2 x_vec ;

		// for each block of the block_row, mult
		for ( int i=rowptr_bounds.x; i<rowptr_bounds.y; i++ ) { 

			row0 = mat0[i] ;
			row1 = mat1[i] ;
			ci = colind[i] ;
			#ifndef USE_TEXTURE
			x_vec = x[ci] ;
			#else
			x_vec = fetch_x(ci,x) ;
			#endif

			res.x += row0.x*x_vec.x+row0.y*x_vec.y ;
			res.y += row1.x*x_vec.x+row1.y*x_vec.y ;
		}
		b[index] = res ;
	}
}


//-------------------------------------------------------------------------------------//
// Main entry point for the CNC plugin solver : it solves Ax=b with                    //
// the matrix A is NL SparseMatrix Format                                              //
// the rhs b is double *                                                               //
// the vector solution x is double *                                                   //
// solvertype must belongs to { FLOAT_CRS | DOUBLE_CRS | FLOAT_BCRS2 | DOUBLE_BCRS2    //
//                              FLOAT_ELL | DOUBLE_ELL | FLOAT_HYB | DOUBLE_HYB }      //
//-------------------------------------------------------------------------------------//

NLuint cnc_solve_cg (
        NLSparseMatrix *A,
        NLdouble *b,
        NLdouble *x,
		NLuint nb_iter_max,
		NLdouble epsilon,
		NLint solver_type ) {
	
	//check cuda device support
	CNC_ASSERT ( CNCCheckDevice (), "No CUDA device found, please check the compatibility of your hardware with CUDA on NVIDIA's website" ) ;
	CNC_ASSERT ( A->n == A->m, "Non-square sparse matrix unsupported" ) ;
	CNC_ASSERT ( nlSolverIsCNC(solver_type), "Wrong solver type" ) ;
	CNC_ASSERT ( nb_iter_max>=1, "Wrong number if iterations" ) ;
	CNC_ASSERT ( epsilon>0.0, "Wrong convergence threshold" ) ;
        CNC_ASSERT ( CNCConfigureDevice (), "CUDA device could not be configured properly");

	NLuint block_size = 1;

    switch(solver_type) {
        case NL_CNC_FLOAT_CRS:
        case NL_CNC_DOUBLE_CRS:
        case NL_CNC_FLOAT_ELL:
        case NL_CNC_DOUBLE_ELL:
        case NL_CNC_FLOAT_HYB:
        case NL_CNC_DOUBLE_HYB: {
            block_size = 1;
        } break ;
        case NL_CNC_FLOAT_BCRS2:
        case NL_CNC_DOUBLE_BCRS2: {
            block_size = 2;
        } break ;
        default:
            nl_assert_not_reached ;
        
    }

    printf ( "############################################################\n" ) ;
	printf ( "Start Init CG Solver: data conversion, allocation and upload\n" ) ;
	printf ( "------------------------------------------------------------\n" ) ;
	printf ( "max iter: %d\ntolerance: %e\nblock_size: %dx%d\n",
			 nb_iter_max, epsilon, block_size, block_size ) ;
	printf ( "size vector: %d\n", A->n ) ;

	NLuint val_ret=static_cast<NLuint>(-1);
    

    switch(solver_type) {
        case NL_CNC_FLOAT_CRS:
        case NL_CNC_FLOAT_BCRS2: {
            printf("using SINGLE precision floating point.\n");
            CNCArray1d<float>  array_x( A->n );
            CNCArray1d<float>  array_b( A->n );	
	
	        std::transform(x,x+A->n,array_x.data(),convert<double,float>());
	        std::transform(b,b+A->n,array_b.data(),convert<double,float>());
	
	        if ( block_size == 1 ) {
                CNCSparseMatrixCRS<float> smcrs ;
		        convert_matrix<float> ( A, smcrs, false ) ;
		        val_ret = solve_cg_internal<CNCSparseMatrixCRS<float>, float> ( 
                    smcrs, array_b, array_x, nb_iter_max, epsilon ) ; 
	        } else if ( block_size == 2 ) {
		        CNCSparseMatrixBCRS<float,2,2> smbcrs2x2 ;
		        convert_matrix<float, 2, 2> ( A, smbcrs2x2 ) ;
		        val_ret = solve_cg_internal<CNCSparseMatrixBCRS<float,2,2> ,float > (  
			        smbcrs2x2, array_b, array_x, nb_iter_max, epsilon ) ; 
	        } else {
		        printf ( "Wrong Block size\n" ) ; // you should never reach this point...
	        }
		
	        std::transform(array_x.data(),array_x.data()+A->n,x,convert<float,double>());
	        std::transform(array_b.data(),array_b.data()+A->n,b,convert<float,double>());
        } break ;
        case NL_CNC_FLOAT_ELL:{
            printf("using SINGLE precision floating point.\n");
            CNCArray1d<float>  array_x( A->n );
            CNCArray1d<float>  array_b( A->n );	
	
	        std::transform(x,x+A->n,array_x.data(),convert<double,float>());
	        std::transform(b,b+A->n,array_b.data(),convert<double,float>());
	
            CNCSparseMatrixELL<float> smell ;
		    convert_matrix<float> ( A, smell) ;
		    val_ret = solve_cg_internal<CNCSparseMatrixELL<float>, float> ( 
		        smell, array_b, array_x,  nb_iter_max, epsilon ) ; 
	        std::transform(array_x.data(),array_x.data()+A->n,x,convert<float,double>());
	        std::transform(array_b.data(),array_b.data()+A->n,b,convert<float,double>());

        } break ;
        case NL_CNC_FLOAT_HYB:{
            printf("using SINGLE precision floating point.\n");
            CNCArray1d<float>  array_x( A->n );
            CNCArray1d<float>  array_b( A->n );	
	
	        std::transform(x,x+A->n,array_x.data(),convert<double,float>());
	        std::transform(b,b+A->n,array_b.data(),convert<double,float>());
	
            CNCSparseMatrixHYB<float> smhyb ;
		    convert_matrix<float> ( A, smhyb) ;
		    val_ret = solve_cg_internal<CNCSparseMatrixHYB<float>, float> ( 
		        smhyb, array_b, array_x,  nb_iter_max, epsilon ) ; 
	        
		    std::transform(array_x.data(),array_x.data()+A->n,x,convert<float,double>());
	        std::transform(array_b.data(),array_b.data()+A->n,b,convert<float,double>());

        } break ;
        case NL_CNC_DOUBLE_ELL: {
            CNC_ASSERT ( CNCCheckDeviceDoubleSupport() , "No support for double precision floating point on the CUDA device found" ) ;
            printf("Nice ! using DOUBLE precision floating point.\n");
            CNCArray1d<double>  array_x( A->n );
            CNCArray1d<double>  array_b( A->n );	
	
	        std::copy(x,x+A->n,array_x.data());
	        std::copy(b,b+A->n,array_b.data());
	
		    CNCSparseMatrixELL<double> smell ;
	    	convert_matrix<double> ( A, smell) ;
		    val_ret = solve_cg_internal<CNCSparseMatrixELL<double>, double> ( 
		        smell, array_b, array_x,  nb_iter_max, epsilon ) ; 
	
	        std::copy(array_x.data(),array_x.data()+A->n,x);
	        std::copy(array_b.data(),array_b.data()+A->n,b);

        } break ;
	    case NL_CNC_DOUBLE_HYB: {
            CNC_ASSERT ( CNCCheckDeviceDoubleSupport() , "No support for double precision floating point on the CUDA device found" ) ;
            printf("Nice ! using DOUBLE precision floating point.\n");
            CNCArray1d<double>  array_x( A->n );
            CNCArray1d<double>  array_b( A->n );	
	
	        std::copy(x,x+A->n,array_x.data());
	        std::copy(b,b+A->n,array_b.data());
	
	        CNCSparseMatrixHYB<double> smhyb ;
		    convert_matrix<double> ( A, smhyb ) ;
		    val_ret = solve_cg_internal<CNCSparseMatrixHYB<double>, double> ( 
	            smhyb, array_b, array_x,  nb_iter_max, epsilon ) ; 
	
	        std::copy(array_x.data(),array_x.data()+A->n,x);
	        std::copy(array_b.data(),array_b.data()+A->n,b);

        } break ;
        case NL_CNC_DOUBLE_CRS:
        case NL_CNC_DOUBLE_BCRS2: {
            CNC_ASSERT ( CNCCheckDeviceDoubleSupport() , "No support for double precision floating point on the CUDA device found" ) ;
        
            printf("Nice ! using DOUBLE precision floating point.\n");
            CNCArray1d<double>  array_x( A->n );
            CNCArray1d<double>  array_b( A->n );	
	
	        std::copy(x,x+A->n,array_x.data());
	        std::copy(b,b+A->n,array_b.data());
	
	        if ( block_size == 1 ) {
		        CNCSparseMatrixCRS<double> smcrs ;
		        convert_matrix<double> ( A, smcrs, false ) ;
		        val_ret = solve_cg_internal<CNCSparseMatrixCRS<double>, double> ( 
			        smcrs, array_b, array_x,  nb_iter_max, epsilon ) ; 
            } else if ( block_size == 2 ) {
		        CNCSparseMatrixBCRS<double,2,2> smbcrs2x2 ;
		        convert_matrix<double, 2, 2> ( A, smbcrs2x2 ) ;
		        val_ret = solve_cg_internal<CNCSparseMatrixBCRS<double,2,2> ,double> ( 
			        smbcrs2x2, array_b, array_x, nb_iter_max, epsilon ) ; 
	        } else {
		        printf ( "Wrong Block size\n" ) ; // you should never reach this point...
	        }
	        std::copy(array_x.data(),array_x.data()+A->n,x);
	        std::copy(array_b.data(),array_b.data()+A->n,b);
        } break ;
        default:
            nl_assert_not_reached ;
    }

	return val_ret;
}


