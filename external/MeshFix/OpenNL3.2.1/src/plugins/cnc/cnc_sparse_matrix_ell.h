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

#ifndef CNC_SPARSE_MATRIX_ELL_H
#define CNC_SPARSE_MATRIX_ELL_H

#include "cnc_utils.h"
#include "cnc_cublas_utils.h"
#include "cnc_arrays.h"
#include "cnc_kernels.h"
#include "cnc_texture.h"
#include <ostream>
#include <iostream>
#include <algorithm>
#include <iterator>

extern "C" {
#include <NL/nl_matrix.h>
}


class CNCSparseMatrixPatternELL {
public:  
    typedef unsigned int index_type;
    typedef CNCArray1d<index_type> IndexArray ;

    CNCSparseMatrixPatternELL() {
        num_rows          = 0 ;
        num_cols          = 0 ;
        num_nonzeros      = 0 ;
        num_cols_per_row  = 0 ;
        stride            = 0 ;
    }

    void clear() {
        num_rows          = 0 ;
        num_cols          = 0 ;
        num_nonzeros      = 0 ;
        num_cols_per_row  = 0 ;
        stride            = 0 ;
        indices.clear() ;
    }


    // number of rows 
    unsigned int m() const { return num_rows ; }

    /** number of columns */
    unsigned int n() const { return num_cols ; }

    /** number of non-zero coefficients */
    unsigned int nnz() const { return num_nonzeros ; }
    
    /** number of num_cols_per_rows */
    unsigned int ncpr() const { return num_cols_per_row ; }


public:
    index_type num_rows;
    index_type num_cols;
    index_type num_nonzeros;
    index_type num_cols_per_row;
    index_type stride;
    
    IndexArray indices;

};

template<typename TypeMatrix> 
class CNCSparseMatrixELL : public CNCSparseMatrixPatternELL {
public:

	typedef TypeMatrix   type_value;
	typedef CNCArray1d<TypeMatrix> CoeffArray ;

    inline CNCSparseMatrixELL()  { }
       
    template <typename TypeVector>
	inline void mult(const TypeVector* x, 
	                 TypeVector* y,
			         const unsigned int vec_size) const ;
        
	template <typename TypeVector>
    inline void mult(const CNCArray1d<TypeVector>& x, 
	                 CNCArray1d<TypeVector>& y) const {
        mult(x.data(), y.data(), x.size()) ;
    }

    inline unsigned int data_size() const {
        return a.size() ;
	}    
        
	inline TypeMatrix * diagonal() {
	    return diag.data() ;
	}    
	
	inline const TypeMatrix * diagonal() const {
	    return diag.data() ;
	}    
        
	inline void gpu_allocate_and_upload () ;
    inline void gpu_deallocate () ;
        
	
	template <typename TypeVector>
	inline void gpu_mult ( const TypeVector * x, 
	                       TypeVector* y, 
			               unsigned int vec_size) ;


    //inline void print(std::ostream& out) ;
	inline void print (); 

    // arrays to store data and diag
	CoeffArray a ;
    CoeffArray diag ;

	// GPU data
	unsigned int  * gpu_indices_ ;
	TypeMatrix * gpu_mat_ ;
} ;

//---------------------------------------------------------------------------//
// printf function for ELL sparse matrix                                     //
//---------------------------------------------------------------------------//


template <typename TypeMatrix>
inline void CNCSparseMatrixELL<TypeMatrix>::print () {
    for(unsigned int i=0; i< m(); i++) {
        for(unsigned int jj=0 ; jj< ncpr() ; jj++) {
            unsigned int index= stride * jj + i;
	        if(a[index] != 0.0){
	            printf ( "%d %d %f\n", i, indices[index], a[index] ) ;
            }
        }
    }
}

//---------------------------------------------------------------------------//
// Matrix-Vector product for ELL sparse matrix on CPU     	    	         //
//---------------------------------------------------------------------------//


template <typename TypeMatrix> 
template <typename TypeVector>
inline void CNCSparseMatrixELL<TypeMatrix>::mult(
        const TypeVector* x, 
        TypeVector* y,
        const unsigned int vec_size) const {
   
    for(unsigned int i=0; i< num_rows; i++) {
        TypeVector sum = 0.0 ;
        for(unsigned int j=0 ; j< num_cols_per_row ; j++) {
            unsigned int index= stride * j + i; 
	        if (indices[index] < vec_size){
	            sum += static_cast<TypeVector>(a[index]) * x[indices[index]] ;
            }
        }
        if (i < vec_size){
            y[i] = sum ;
        }
   }
}

//---------------------------------------------------------------------------//
// Allocate and upload a ELL sparse matrix on GPU memory					 //
//---------------------------------------------------------------------------//

template <typename TypeMatrix>
inline void CNCSparseMatrixELL<TypeMatrix>::gpu_allocate_and_upload () {

	// allocate and upload non-zero matrix coefficients
	new_device_vector ( &gpu_mat_, a.size()+16 ) ;
	copy_host_to_device_vector( a, gpu_mat_ ) ;

	// allocate and upload indices
	new_device_vector ( &gpu_indices_, indices.size()+16 ) ;
	copy_host_to_device_vector ( indices , gpu_indices_ ) ;

}

//---------------------------------------------------------------------------//
// Free graphics memory used to store a ELL sparse matrix					 //
//---------------------------------------------------------------------------//

template <typename TypeMatrix> 
inline void CNCSparseMatrixELL<TypeMatrix>::gpu_deallocate () {

	delete_device_vector ( &gpu_mat_ ) ;
	delete_device_vector ( &gpu_indices_ ) ;
}


//---------------------------------------------------------------------------//
// GPU ELL Matrix / Vector Multiply											 //
//---------------------------------------------------------------------------//
template <typename TypeMatrix>
template <typename TypeVector> 
inline void CNCSparseMatrixELL<TypeMatrix>::gpu_mult (
		const TypeVector * x, 
		TypeVector * y,
		unsigned int vec_size) {
    
    #ifdef USE_TEXTURE
    // store the x vector in texture
    bind_x(x);
    #endif
    
    CNC_SIMPLE_ASSERT(
        (num_cols == num_rows ) && (num_cols == vec_size),
        "matrix not square or vec_size != num_cols!");
    
    const unsigned int BLOCK_SIZE = 256;
    const dim3 grid = make_large_grid(num_rows, BLOCK_SIZE);
    
    CNCMatELL_Nathan_Bell<TypeMatrix> <<<grid, BLOCK_SIZE>>> (
        num_rows,
        num_cols, 
        num_cols_per_row, 
        stride, 
        gpu_mat_, 
        gpu_indices_, 
        x, 
        y); 
    
    #ifdef USE_TEXTURE
    // store the x vector in texture
    unbind_x(x);
    #endif
}


//---------------------------------------------------------------------------//
//  Conversion NLSparseMatrix to CNCSparseMatrixELL<TypeMatrix> : definition //	
//---------------------------------------------------------------------------//

template <typename TypeMatrix> 
inline void convert_matrix( NLSparseMatrix * rhs, 
                            CNCSparseMatrixELL<TypeMatrix>& A, 
                            const unsigned int alignment = 32 ) {

	CNC_SIMPLE_ASSERT ( rhs->n == rhs->m, "In CNC, matrices must be square!" ) ;
    A.num_cols = rhs->n ;
    A.num_rows = rhs->m ;
	

	unsigned int num_nonzeros = nlSparseMatrixNNZ(rhs) ;
        
    // seeking for the max number of nnz elements per row
    unsigned int max_nnz_per_rows = 0;
    for(unsigned int i=0; i<rhs->m; i++) {
        max_nnz_per_rows = std::max(max_nnz_per_rows,rhs->row[i].size); 
    }

    // compute the stride 
    A.stride = alignment * ((A.num_rows + alignment - 1)/ alignment);
       
    A.num_cols_per_row = max_nnz_per_rows;
	unsigned int common_size = A.stride * A.num_cols_per_row ;
	
	// allocate 
	A.indices.allocate( common_size) ;
	A.a.allocate( common_size) ;
    A.diag.allocate( A.num_rows ) ;
       
    //  filling with 0 
	std::fill( A.indices.data(), A.indices.data() + common_size, 0);
    std::fill(       A.a.data(),       A.a.data() + common_size, static_cast<TypeMatrix>(0));
	std::fill(    A.diag.data(),    A.diag.data() + A.num_rows , static_cast<TypeMatrix>(0));
	
	for(unsigned int i=0; i<rhs->m; i++) {
	    
	    NLRowColumn * R = &(rhs->row[i]) ;
	    for(unsigned int jj=0; jj<R->size; jj++) {
                 
		    // index of (i,jj) in A
		    unsigned int index = A.stride* jj + i;
		 
		    A.a[index]       = static_cast<TypeMatrix>(R->coeff[jj].value) ;
            A.indices[index] =                         R->coeff[jj].index  ;
            // if diagonal element we store it in diag
	        if (i == R->coeff[jj].index){ 
		        A.diag[i] = static_cast<TypeMatrix>(R->coeff[jj].value);
	        }
        }
    }
    std::cout << "average nnz per rows = "  << num_nonzeros / rhs->n << std::endl
              << "max     nnz per rows = "  << max_nnz_per_rows << std::endl;
}

#endif
