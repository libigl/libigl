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

#ifndef CNC_SPARSE_MATRIX_CRS_H
#define CNC_SPARSE_MATRIX_CRS_H

#include "cnc_utils.h"
#include "cnc_cublas_utils.h"
#include "cnc_arrays.h"
#include "cnc_kernels.h"
#include <ostream>
#include <iostream>
#include <algorithm>

extern "C" {
#include <NL/nl_matrix.h>
}



// -------------------------------------------------------------------------- //
// CRS Sparse matrix data structure with CPU and GPU storage				  //
// -------------------------------------------------------------------------- //

class CNCSparseMatrixPatternCRS {
public:
    typedef CNCArray1d<unsigned int> IndexArray ;
    enum { NO_INDEX = ~0 } ;

    CNCSparseMatrixPatternCRS() {
        N = 0 ;
        symmetric_storage = false ;
    }

    void clear() {
        N = 0 ;
        symmetric_storage = false ;
        colind.clear() ;
        rowptr.clear() ;
    }


    /** number of rows */
    unsigned int m() const { return rowptr.size() - 1 ; }

    /** number of columns */
    unsigned int n() const { return N ; }

    /** number of non-zero coefficients */
    unsigned int nnz() const { return colind.size() ; }

    /** number of non-zero coefficients in row i */
    unsigned int row_nnz(unsigned int i) const { 
        return rowptr[i+1] - rowptr[i] ;
    }

public:
    IndexArray colind ;
    IndexArray rowptr ;

    unsigned int N ;
    bool symmetric_storage ;

} ;

//---------------------------------------------------------------------------//

template<typename TypeMatrix> 
class CNCSparseMatrixCRS : public CNCSparseMatrixPatternCRS {
public:

	typedef TypeMatrix type_value;
	typedef CNCArray1d<TypeMatrix> CoeffArray ;

    inline CNCSparseMatrixCRS() : separate_diag(true) { }
        
	template <typename TypeVector>
    inline void mult(const TypeVector* x, 
	                 TypeVector* y,
			         const unsigned int vec_size) const ;
        
	template <typename TypeVector>
    inline void mult(const CNCArray1d<TypeVector>& x, 
	                 CNCArray1d<TypeVector>& y) const {
        mult( x.data(), y.data(), x.size()) ;
    }
    inline void print(std::ostream& out) ;


    inline unsigned int data_size() const {
	    return a.size();
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
	inline void gpu_mult ( const TypeVector * x, TypeVector * y, 
	                       unsigned int vec_size) ;


	inline void print () {
	    for(unsigned int i=0; i<m(); i++) {
		    for(unsigned int jj=rowptr[i]; jj<rowptr[i+1]; jj++) {
			    if(a[jj] != 0.0) {
				    unsigned int j = colind[jj] ;
				    printf ( "%d %d %f\n", i, j, a[jj] ) ;
				    if(symmetric_storage && j != i) {
					    printf ( "%d %d %f\n", j, i, a[jj] ) ;
				    }
			    }
		    }
	    }
	    if(separate_diag) {
		    for(unsigned int i=0; i<diag.size(); i++) {
			    if(diag[i] != 0.0) {
				    printf ( "%d %d %f\n", i, i, diag[i] ) ;
			    }
		    }
	    }
    }

    // CPU data
    CoeffArray a ;
    CoeffArray diag ;
    bool separate_diag ;

	// GPU data
	uint2 * gpu_redundant_rp_ ;
	unsigned int * gpu_ci_ ;
	TypeMatrix * gpu_mat_ ;
} ;

//---------------------------------------------------------------------------//
// CRS Sparse Matrix-Vector product on CPU                                   //
//---------------------------------------------------------------------------//

template <typename TypeMatrix> 
template <typename TypeVector> 
inline void CNCSparseMatrixCRS<TypeMatrix>::mult(
        const  TypeVector * x, 
        TypeVector* y,
		const unsigned int vec_size) const {

    if(separate_diag) {
        if(symmetric_storage) {
            // Symmetric storage and separate diagonal
            for(unsigned int i=0; i<diag.size(); i++) {
                y[i] = x[i] * static_cast<TypeVector>(diag[i]) ;
            }
            for(unsigned int i=diag.size(); i<m(); i++) {
                y[i] = 0.0 ;
            }
            for(unsigned int i=0; i<m(); i++) {
                for(unsigned int j=rowptr[i] ; j<rowptr[i+1]; j++) {
                    if ((colind[j] < vec_size) & (i < vec_size)){
		                y[i] += static_cast<TypeVector>(a[j]) * x[colind[j]] ;
                        y[colind[j]] += static_cast<TypeVector>(a[j]) * x[i] ;                        
                    }
	            }
            }
        } else {
            // No symmetric storage and separate diagonal
            for(unsigned int i=0; i<m(); i++) {
                TypeVector sum = 0.0 ;
                for(unsigned int j=rowptr[i] ; j<rowptr[i+1]; j++) {
                    sum += static_cast<TypeVector>(a[j]) * x[colind[j]] ;
                }
                y[i] = sum ;
            }
            for(unsigned int i=0; i<diag.size(); i++) {
                y[i] += x[i] * static_cast<TypeVector>(diag[i]) ;
            }
        }
    } else {
        if(symmetric_storage) {
            // Symmetric storage, no separate diagonal
	        memset ( y, 0, sizeof(TypeVector) ) ;
            for(unsigned int i=0; i<m(); i++) {
                for(unsigned int j=rowptr[i] ; j<rowptr[i+1]; j++) {
                    y[i] += static_cast<TypeVector>(a[j]) * x[colind[j]] ;
                    if(colind[j] != i) {
                        y[colind[j]] += static_cast<TypeVector>(a[j]) * x[i] ;                        
                    }
                }
            }
        } else {
            // No symmetric storage, no separate diagonal
            for(unsigned int i=0; i<m(); i++) {
                TypeVector sum = 0.0 ;
                for(unsigned int j=rowptr[i] ; j<rowptr[i+1]; j++) {
                    sum += static_cast<TypeVector>(a[j]) * x[colind[j]] ;
                }
                y[i] = sum ;
            }
        }
    }
}


//---------------------------------------------------------------------------//
// Allocate and upload a CRS sparse matrix on GPU memory                     //
//---------------------------------------------------------------------------//

template <typename TypeMatrix> 
inline void CNCSparseMatrixCRS<TypeMatrix>::gpu_allocate_and_upload () {

	// allocate and upload non-zero matrix coefficients
	new_device_vector( &gpu_mat_, a.size()+16 );
	copy_host_to_device_vector(  a, gpu_mat_ );

	// allocate and upload rowptr and colind
	new_device_vector( &gpu_ci_, colind.size()+16 ) ;
	new_device_vector( &gpu_redundant_rp_, (rowptr.size()-1)+16 ) ;

	uint2 * cpu_redundant_rp ;
	new_host_vector( &cpu_redundant_rp, rowptr.size()-1 +16) ;

	for (unsigned int i=0; i<rowptr.size()-1; i++) {
		cpu_redundant_rp[i].x = rowptr.data()[i  ] ;
		cpu_redundant_rp[i].y = rowptr.data()[i+1] ;
	}

	copy_host_to_device_vector (
	    cpu_redundant_rp,
	    gpu_redundant_rp_,
	    rowptr.size()-1+16 ) ;
	copy_host_to_device_vector ( colind , gpu_ci_) ;

	delete_host_vector ( &cpu_redundant_rp ) ;
}


//---------------------------------------------------------------------------//
// Free graphics memory used to store a CRS sparse matrix                    //
//---------------------------------------------------------------------------//

template <typename TypeMatrix> 
inline void CNCSparseMatrixCRS<TypeMatrix>::gpu_deallocate () {
	delete_device_vector ( &gpu_mat_ ) ;
	delete_device_vector ( &gpu_ci_ ) ;
	delete_device_vector ( &gpu_redundant_rp_ ) ;
}


//---------------------------------------------------------------------------//
// GPU CRS Matrix / Vector Multiply                                          //
//---------------------------------------------------------------------------//

template <typename TypeMatrix> template <typename TypeVector> 
inline void CNCSparseMatrixCRS<TypeMatrix>::gpu_mult ( const TypeVector * x, 
                                                       TypeVector * y,
		                                               unsigned int vec_size) {
    #ifdef USE_TEXTURE
    // store the x vector in texture
    bind_x(x);
    #endif
    
    #ifndef USE_SHARED_CRS 
    const unsigned int BLOCK_SIZE = 256;
    const dim3 grid = make_large_grid(vec_size, BLOCK_SIZE);
    CNCMat1x1VecMultKernel<TypeMatrix,TypeVector> <<<grid, BLOCK_SIZE>>>( 
         gpu_mat_,
         a.size(),
         gpu_redundant_rp_,
         rowptr.size(),
         gpu_ci_,
         colind.size(),
         x,
         y,
         vec_size); 
    #else
    const unsigned int NB_BLOCK_SIZE = 128;
    const unsigned int WARPS_PER_BLOCK = NB_BLOCK_SIZE / WARP_SIZE;
    const unsigned int MAX_BLOCKS = MAX_THREADS / NB_BLOCK_SIZE;
    const unsigned int NUM_BLOCKS = std::min(MAX_BLOCKS, DIVIDE_INTO(m(), WARPS_PER_BLOCK));
    CNCCRSNathan_Bell<TypeMatrix,TypeVector,NB_BLOCK_SIZE> <<<NUM_BLOCKS, NB_BLOCK_SIZE>>> (
        gpu_mat_,
        a.size(),
        gpu_redundant_rp_,
        rowptr.size(),
        gpu_ci_,
        colind.size(),
        x,
        y,
        vec_size); 
    #endif

    #ifdef USE_TEXTURE
    // store the x vector in texture
    unbind_x(x);
    #endif
}

//---------------------------------------------------------------------------//
//  Conversion NLSparseMatrix to CNCSparseMatrixCRS<T> : definition          //	
//---------------------------------------------------------------------------//

template <typename TypeMatrix> 
inline void convert_matrix( NLSparseMatrix * rhs, 
                            CNCSparseMatrixCRS<TypeMatrix>& A, 
                            bool separate_diag ) {
    A.separate_diag = separate_diag ;

    A.symmetric_storage = static_cast<bool>(rhs->storage &  NL_MATRIX_STORE_SYMMETRIC);
    A.N = rhs->n ;
    A.rowptr.allocate(rhs->m + 1) ;
    unsigned int nnz = nlSparseMatrixNNZ(rhs) ;
        
	if(separate_diag) {
        nnz -= rhs->diag_size;
    }
        
	// allocation 
	A.colind.allocate(nnz) ;
    A.a.allocate(nnz) ;
    A.diag.allocate(rhs->diag_size) ;
        
    unsigned int max_nnz_per_rows = 0;
	unsigned int cur = 0 ;
    for(unsigned int i=0; i<rhs->m; i++) {
        A.rowptr[i] = cur ;
        NLRowColumn * R = &(rhs->row[i]) ;
        max_nnz_per_rows = std::max(max_nnz_per_rows,R->size); 
	    for(unsigned int jj=0; jj<R->size; jj++) {
	     	if(!separate_diag || (R->coeff[jj].index != i )) {
                A.a[cur] = static_cast<TypeMatrix>(R->coeff[jj].value) ;
                A.colind[cur] = R->coeff[jj].index ;
                cur++ ;
            }
        }
    }
        
	A.rowptr[rhs->m] = nnz ;
        
	// we retrieve the diagonal
	for(unsigned int i=0; i<rhs->diag_size; i++) {
        A.diag[i] = static_cast<TypeMatrix>(rhs->diag[i]) ;
    }
        
	std::cout << "average nnz per rows = "  << nnz / rhs->n << std::endl
	          << "max     nnz per rows = "  << max_nnz_per_rows << std::endl;
}

#endif
