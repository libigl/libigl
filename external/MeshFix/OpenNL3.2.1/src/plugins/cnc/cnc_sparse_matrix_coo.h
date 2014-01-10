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


#ifndef CNC_SPARSE_MATRIX_COO_H
#define CNC_SPARSE_MATRIX_COO_H

#include "cnc_utils.h"
#include "cnc_cublas_utils.h"
#include "cnc_arrays.h"
#include "cnc_kernels.h"
#include "cnc_texture.h"
#include "cnc_sparse_matrix_crs.h"
#include <ostream>
#include <iostream>
#include <algorithm>
#include <iterator>

extern "C" {
#include <NL/nl_matrix.h>
}


class CNCSparseMatrixPatternCOO {
public:  
    typedef unsigned int index_type;
    typedef CNCArray1d<index_type> IndexArray ;

    CNCSparseMatrixPatternCOO() {
        num_rows          = 0 ;
        num_cols          = 0 ;
        num_nonzeros      = 0 ;
    }

    void clear() {
        num_rows          = 0 ;
        num_cols          = 0 ;
        num_nonzeros      = 0 ;
        I.clear() ;
        J.clear() ;
    }

    // number of rows 
    unsigned int m() const { return num_rows ; }

    /** number of columns */
    unsigned int n() const { return num_cols ; }

    /** number of non-zero coefficients */
    unsigned int nnz() const { return num_nonzeros ; }
    

public:
    index_type num_rows;
    index_type num_cols;
    index_type num_nonzeros;
    
    IndexArray I;
    IndexArray J;

};

template<typename TypeMatrix> 
class CNCSparseMatrixCOO : public CNCSparseMatrixPatternCOO {
public:

	typedef TypeMatrix   type_value;
	typedef CNCArray1d<TypeMatrix> CoeffArray ;

    inline CNCSparseMatrixCOO()  { }
   
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
	    return V.size() ;
	}    

    inline TypeMatrix * diagonal() {
        return diag.data();
	}

    inline const TypeMatrix * diagonal() const {
        return diag.data();
	}
        
	inline void gpu_allocate_and_upload () ;
    inline void gpu_deallocate () ;
    
	template <typename TypeVector>
	inline void gpu_mult ( const TypeVector * x, 
	                       TypeVector* y, 
			               unsigned int vec_size) ;

	inline void print (); 

    // arrays to store data and diag
	CoeffArray V ;
    CoeffArray diag ;

	// GPU data
	unsigned int * gpu_I_ ;
	unsigned int * gpu_J_ ;
	TypeMatrix * gpu_V_ ;
};

//---------------------------------------------------------------------------//
// printf function for COO sparse matrix                                     //
//---------------------------------------------------------------------------//


template <typename TypeMatrix> 
inline void CNCSparseMatrixCOO<TypeMatrix>::print () {
    for(unsigned int i=0; i< nnz() ; i++){ 
	    printf ( "%d %d %f\n", I[i], J[i], V[i] ) ;
    }
}

//---------------------------------------------------------------------------//
// Matrix-Vector product for COO sparse matrix on CPU                        //
//---------------------------------------------------------------------------//


template <typename TypeMatrix> 
template <typename TypeVector>
inline void CNCSparseMatrixCOO<TypeMatrix>::mult(
        const TypeVector* x, 
        TypeVector* y,
		const unsigned int vec_size) const {
    std::fill(y,y + vec_size, static_cast<TypeMatrix>(0));
    for(unsigned int i=0; i< num_nonzeros ; i++) {
        if (I[i] < vec_size){
	        y[I[i]] += static_cast<TypeVector>(V[i]) * x[J[i]] ;
        }
    }
}

//---------------------------------------------------------------------------//
// Allocate and upload a COO sparse matrix on GPU memory					 //
//---------------------------------------------------------------------------//

template <typename TypeMatrix>
inline void CNCSparseMatrixCOO<TypeMatrix>::gpu_allocate_and_upload () {

	// allocate the I,J,V triplet 
	new_device_vector ( &gpu_V_, V.size()+16 ) ;
	new_device_vector ( &gpu_I_, I.size()+16 ) ;
	new_device_vector ( &gpu_J_, J.size()+16 ) ;
    // and upload!
	copy_host_to_device_vector ( V , gpu_V_ ) ;
	copy_host_to_device_vector ( I , gpu_I_ ) ;
	copy_host_to_device_vector ( J , gpu_J_ ) ;

}

//---------------------------------------------------------------------------//
// Free graphics memory used to store a COO sparse matrix                    //
//---------------------------------------------------------------------------//

template <typename TypeMatrix> 
inline void CNCSparseMatrixCOO<TypeMatrix>::gpu_deallocate () {

	delete_device_vector ( &gpu_V_ ) ;
	delete_device_vector ( &gpu_I_ ) ;
	delete_device_vector ( &gpu_J_ ) ;
}


//---------------------------------------------------------------------------//
// GPU COO Matrix / Vector Multiply : Taken Directly From Nathan Bell spmv   //
//---------------------------------------------------------------------------//
template <typename TypeMatrix>
template <typename TypeVector> 
inline void CNCSparseMatrixCOO<TypeMatrix>::gpu_mult (
		const TypeVector * x, 
		TypeVector * y,
		unsigned int vec_size) {
 
    // fill y with zeroes 
    const unsigned int TMP_BLOCK_SIZE = 256;
    const dim3 grid = make_large_grid(vec_size, TMP_BLOCK_SIZE);
    CNCfill <<<grid,TMP_BLOCK_SIZE>>>(y,vec_size,static_cast<TypeVector>(0));
    
    if(num_nonzeros == 0){
        // empty matrix
        return;
    } else if (num_nonzeros < WARP_SIZE){
        // small matrix
        spmv_coo_serial_kernel<TypeMatrix,TypeVector> <<<1,1>>> (
            num_nonzeros,
            gpu_I_,
            gpu_J_,
            gpu_V_,
            x,
            y);
        return;
    }

    
    const unsigned int BLOCK_SIZE      = 256;
    const unsigned int MAX_BLOCKS      = MAX_THREADS / (2 * BLOCK_SIZE);
    const unsigned int WARPS_PER_BLOCK = BLOCK_SIZE / WARP_SIZE;

    const unsigned int num_units  = num_nonzeros / WARP_SIZE; 
    const unsigned int num_warps  = std::min(
        num_units,
        WARPS_PER_BLOCK * MAX_BLOCKS);
    const unsigned int num_blocks = DIVIDE_INTO(num_warps, WARPS_PER_BLOCK);
    const unsigned int num_iters  = DIVIDE_INTO(num_units, num_warps);
    
    const unsigned int interval_size = WARP_SIZE * num_iters;

    // do the last few nonzeros separately (fewer than WARP_SIZE elements)
    const unsigned int tail = num_units * WARP_SIZE; 

    const unsigned int active_warps = (interval_size == 0) ? 0 : DIVIDE_INTO(tail, interval_size);

    #ifdef USE_TEXTURE
    bind_x(x);
    #endif

    unsigned int * temp_rows ; 
    TypeVector * temp_vals ; 

    new_device_vector( &temp_rows, active_warps );
    new_device_vector( &temp_vals, active_warps );

    spmv_coo_flat_kernel<TypeMatrix,TypeVector, BLOCK_SIZE> <<<num_blocks, BLOCK_SIZE>>> 
          ( tail, interval_size, gpu_I_, gpu_J_, gpu_V_, x, y, temp_rows, temp_vals); 
    
    spmv_coo_serial_kernel<TypeMatrix,TypeVector> <<<1,1>>> 
          ( num_nonzeros - tail, gpu_I_ + tail, gpu_J_ + tail, gpu_V_ + tail, x, y);

    spmv_coo_reduce_update_kernel<TypeVector, 512> <<<1, 512>>> 
          ( active_warps, temp_rows, temp_vals, y); 
    delete_device_vector ( &temp_rows );
    delete_device_vector ( &temp_vals );

    #ifdef USE_TEXTURE 
    unbind_x(x);
    #endif
}


//------------------------------------------------------------------------------//
//  Conversion CNCSparseCRSMatrix<TypeMatrix> to CNCSparseMatrixCOO<TypeMatrix> /	
//------------------------------------------------------------------------------//

template <typename TypeMatrix> 
inline void convert_matrix( CNCSparseMatrixCRS<TypeMatrix> &rhs, 
                            CNCSparseMatrixCOO<TypeMatrix> & A ) {

	CNC_SIMPLE_ASSERT ( rhs.n() == rhs.m(), "In CNC, matrices must be square!" ) ;
        
	A.num_nonzeros = rhs.a.size();
	A.num_cols     = rhs.n() ;
    A.num_rows     = rhs.m() ;
	

	// allocation of the triplet I,J,V and diagonal
	A.I.allocate( A.num_nonzeros) ;
	A.J.allocate( A.num_nonzeros) ;
	A.V.allocate( A.num_nonzeros) ;
    A.diag.allocate( A.num_rows ) ;
       
    //  copying data from CRS format to COO format 
    std::copy( rhs.a.data(), rhs.a.data() + A.num_nonzeros, A.V.data() );
    std::copy( rhs.colind.data(), rhs.colind.data() + A.num_nonzeros, A.J.data() );
	std::copy( rhs.diag.data(), rhs.diag.data() + rhs.m(), A.diag.data() );
	
    unsigned int *begin = A.I.data();	
    unsigned int *end;
	for(unsigned int i=0; i< rhs.m() ;  i++) {
        end = begin+ rhs.rowptr[i+1]-rhs.rowptr[i];
        // we "decompress" the crs format in the coo format!
		std::fill( begin, end , i);
        begin = end;
	}
}

#endif
