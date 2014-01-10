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

#ifndef CNC_SPARSE_MATRIX_BCRS_H
#define CNC_SPARSE_MATRIX_BCRS_H

#include "cnc_utils.h"
#include "cnc_cublas_utils.h"
#include "cnc_arrays.h"
#include "cnc_kernels.h"
#include "cnc_texture.h"

#include <algorithm>
#include <set>



// -------------------------------------------------------------------------- //
// BCRS Sparse matrix data structure with CPU and GPU storage	      	      //
// -------------------------------------------------------------------------- //

template <class TypeMatrix, int BLOCM, int BLOCN> 
class CNCGenericSparseMatrixBCRS {
public:
    typedef CNCArray1d<unsigned int> IndexArray ;
    typedef CNCArray1d<TypeMatrix> CoeffArray ;
    enum { BM=BLOCM } ;
    enum { BN=BLOCN } ;
    enum { BLOC_SIZE=BM*BN } ;

    inline CNCGenericSparseMatrixBCRS() : N_(0) { } 

    inline float filling_ratio() const ;

    /** number of bloc rows */
    unsigned int M() const { return rowptr.size() - 1 ; }
    /** number of bloc columns */
    unsigned int N() const { return N_ ; }
    /** number of rows */
    unsigned int m() const { return M() * BM ; } 
    /** number of columns */
    unsigned int n() const { return N() * BN ; } 

    unsigned int number_of_blocs_in_row(unsigned int i) {
        return rowptr[i+1] - rowptr[i] ;
    }

    unsigned int max_number_of_blocs_in_row() {
        unsigned int result = 0 ;
        for(unsigned int i=0; i<M(); i++) {
	        result = (result>number_of_blocs_in_row(i))?(result):(number_of_blocs_in_row(i)) ;
        }
        return result ;
    }
        
    template <typename TypeVector>
    void mult(const TypeVector* x, 
	          TypeVector* y,
              const unsigned int vec_size) const ;

	template <typename TypeVector>
    inline void mult(const CNCArray1d<TypeVector>& x, 
	                 CNCArray1d<TypeVector>& y) const {
        mult( x.data(), y.data(), x.size()) ;
    }

        
	inline unsigned int data_size() const {
	    return a.size() ;
	}    
	
	inline TypeMatrix * diagonal()  {
	    return diag.data() ;
	}    
	
	inline const TypeMatrix * diagonal() const {
	    return diag.data() ;
	}    

    inline void print() const {
        for(unsigned int I=0; I<M(); I++) {
            for(unsigned int JJ=rowptr[I]; JJ<rowptr[I+1]; JJ++) {
                unsigned int ofs = 0 ;
         	    unsigned int J = colind[JJ] ;
           	    for(unsigned int di=0; di<BM; di++) {
     	            for(unsigned int dj=0; dj<BN; dj++) {
     	 	            if (a[JJ * BLOC_SIZE + ofs] != 0.0) {
     		                printf ( "%d %d %f\n",
     		                    (I * BM + di),
     		                    (J * BN + dj),
     		                    a[JJ * BLOC_SIZE + ofs] ) ;
     		            }
     		        ofs++ ;
     	            }
	     	    }
	        }
	    }
	}

    IndexArray rowptr ;
    IndexArray colind ;
    CoeffArray a ;
    CoeffArray diag ;
    unsigned int N_ ;
} ;


//---------------------------------------------------------------------------//
// BLOC CRS Sparse Matric class                                              //
//---------------------------------------------------------------------------//
template <class TypeMatrix, int BLOCM, int BLOCN> 
class CNCSparseMatrixBCRS : 
    public CNCGenericSparseMatrixBCRS<TypeMatrix, BLOCM, BLOCN> {
public:
    typedef TypeMatrix type_value;
    typedef CNCGenericSparseMatrixBCRS<TypeMatrix, BLOCM, BLOCN> superclass ;
    typedef CNCArray1d<unsigned int> IndexArray ;
	
	CNCSparseMatrixBCRS() : 
		gpu_redundant_rp_(NULL),
		gpu_ci_(NULL),
		gpu_mat_(NULL),
		gpu_mat0_(NULL),
		gpu_mat1_(NULL),
		gpu_mat2_(NULL),
		gpu_mat3_(NULL) { }

	inline void gpu_allocate_and_upload () ;
	inline void gpu_deallocate () ;
	
	template <typename TypeVector>
	inline void gpu_mult ( const TypeVector * x, 
	                       TypeVector * y, 
                           unsigned int vec_size) ;
	
    // GPU data
    uint2 * gpu_redundant_rp_ ;
    unsigned int * gpu_ci_ ;
    TypeMatrix * gpu_mat_ ;
    TypeMatrix * gpu_mat0_;
    TypeMatrix * gpu_mat1_;
    TypeMatrix * gpu_mat2_;
    TypeMatrix * gpu_mat3_;
} ;

//---------------------------------------------------------------------------//
// Matrix-Vector product on CPU fror Generic BCRS Sparse Matrix              //
//---------------------------------------------------------------------------//

template <class TypeMatrix, int BM, int BN> 
template <typename TypeVector>
void CNCGenericSparseMatrixBCRS<TypeMatrix,BM,BN>::mult(
        const TypeVector* x, 
        TypeVector* y,
        const unsigned int vec_size) const {
    unsigned int y_base = 0 ;
        
    for(unsigned int I=0; I<M(); I++) {
        TypeVector sum[BM] ;
        for(unsigned int i=0; i<BM; i++) {
            sum[i] = 0 ;
        }
        for(unsigned int JJ=rowptr[I]; JJ<rowptr[I+1]; JJ++) {
            unsigned int ofs = 0 ;
            unsigned int J = colind[JJ] ;
            unsigned int a_base = JJ*BLOC_SIZE ;
            unsigned int x_base = J * BN ;
            for(unsigned int di=0; di<BM; di++) {
                for(unsigned int dj=0; dj<BN; dj++) {
                    if (x_base+dj < vec_size){
                        sum[di] += static_cast<TypeVector>(a[a_base+ofs]) * x[x_base+dj] ;
                    }
                    ofs++ ;
                }
            }
        }
        for(unsigned int i=0; i<BM; i++) {
            if (y_base + i < vec_size ){
                y[y_base + i] = sum[i] ;
            }
        }

        y_base += BM ;
    }   
}


//----------------------------------------------------------------------------//
// Allocation of graphics memory and upload of a 2x2 BCRS float sparse matrix //
//----------------------------------------------------------------------------//

template <> 
inline void CNCSparseMatrixBCRS<float,2,2>::gpu_allocate_and_upload () {

	unsigned int size_matrix = a.size () ;
	
	// allocate and upload non-zero matrix coefficients	
	new_device_vector ( &gpu_mat_, size_matrix  +16 ) ;
	copy_host_to_device_vector ( a , gpu_mat_ ) ;

	// allocate and upload rowptr and colind
	new_device_vector ( &gpu_ci_, colind.size()+16) ;
	new_device_vector ( &gpu_redundant_rp_, (rowptr.size()-1)+16 ) ;

	uint2 * cpu_redundant_rp ;
	new_host_vector( &cpu_redundant_rp, rowptr.size()-1 + 16) ;
	
	for (unsigned int i=0; i<rowptr.size()-1; i++) {
		cpu_redundant_rp[i].x = rowptr.data()[i  ] ;
		cpu_redundant_rp[i].y = rowptr.data()[i+1] ;
	}

	copy_host_to_device_vector ( cpu_redundant_rp, gpu_redundant_rp_, rowptr.size()-1 +16) ;
	copy_host_to_device_vector ( colind , gpu_ci_ ) ;

	delete_host_vector ( &cpu_redundant_rp ) ;
}

//-----------------------------------------------------------------------------//
// Allocation of graphics memory and upload of a 2x2 BCRS double sparse matrix //
//-----------------------------------------------------------------------------//

template <> inline void CNCSparseMatrixBCRS<double,2,2>::gpu_allocate_and_upload () {

	unsigned int size_matrix = a.size () ;
	
	// allocate and upload non-zero matrix coefficients	
	new_device_vector ( &gpu_mat0_, size_matrix/2+16 ) ;
	new_device_vector ( &gpu_mat1_, size_matrix/2+16 ) ;
	
	double * m0 ;
	double * m1 ;
	new_host_vector( &m0, size_matrix/2 + 16 );
	new_host_vector( &m1, size_matrix/2 + 16 );


	// reorder sparse matrix blocs(double)  to fit into two arrays (double2)
	for(unsigned int i=0; i<(unsigned int)(size_matrix/4); i++) {
		m0[i*2  ] = a.data()[i*4   ] ;
		m0[i*2+1] = a.data()[i*4+ 1] ;
		m1[i*2  ] = a.data()[i*4+ 2] ;
		m1[i*2+1] = a.data()[i*4+ 3] ;
	}

	copy_host_to_device_vector ( m0, gpu_mat0_, size_matrix/2+16 ) ;
	copy_host_to_device_vector ( m1, gpu_mat1_, size_matrix/2+16 ) ;

	delete_host_vector( &m0 );
	delete_host_vector( &m1 );
        
	// allocate and upload rowptr on GPU 
	new_device_vector ( &gpu_redundant_rp_, (rowptr.size()-1)+16 ) ;
	new_device_vector ( &gpu_ci_ , colind.size()+16 ) ;
	
	// beware : rowptr is unsigned int but gpu_redundant will be uint2 
	// the point is the  following : (rowptr[i],rowptr[i+1]) are the bounds for the j-index in spVM product    
	// to facilitate the product on GPU we duplicate rowptr[i+1] on CPU
	// in gpu_redundant_rp[i].y and gpu_redundant_rp[i+1].x
	uint2 * cpu_redundant_rp ;
	new_host_vector( &cpu_redundant_rp, rowptr.size()-1 + 16) ;
	
	for (unsigned int i=0; i<rowptr.size()-1; i++) {
		cpu_redundant_rp[i].x = rowptr.data()[i  ] ;
		cpu_redundant_rp[i].y = rowptr.data()[i+1] ;
	}
	
	unsigned int * cpu_ci;
	new_host_vector( &cpu_ci, colind.size() + 16) ;
       
   // copy from colind to cpu_ci
    std::copy(colind.data(), colind.data() + colind.size() , cpu_ci);
        

	copy_host_to_device_vector ( cpu_redundant_rp, gpu_redundant_rp_, rowptr.size()-1 + 16) ;
	copy_host_to_device_vector ( cpu_ci , gpu_ci_ ,colind.size()+16 ) ;

	delete_host_vector( &cpu_ci ) ;
	delete_host_vector( &cpu_redundant_rp ) ;
}


//---------------------------------------------------------------------------//
// Free the gpu memory used to store a BCRS sparse matrix					 //
//---------------------------------------------------------------------------//

template <> inline void CNCSparseMatrixBCRS<float,2,2>::gpu_deallocate () {

	delete_device_vector ( &gpu_ci_ ) ;
	delete_device_vector ( &gpu_redundant_rp_ ) ;
	delete_device_vector ( &gpu_mat_ ) ;
}

//---------------------------------------------------------------------------//
// Free the gpu memory used to store a BCRS double 2x2 matrix	             //
//---------------------------------------------------------------------------//

template <> inline void CNCSparseMatrixBCRS<double,2,2>::gpu_deallocate () {

	delete_device_vector ( &gpu_ci_ ) ;
	delete_device_vector ( &gpu_redundant_rp_ ) ;
	delete_device_vector ( &gpu_mat0_ ) ;
	delete_device_vector ( &gpu_mat1_ ) ;
}


//---------------------------------------------------------------------------//
// GPU BCRS float 2x2 Matrix / Vector Multiply on GPU                        //
//---------------------------------------------------------------------------//

template <> template <typename TypeVector> 
inline void CNCSparseMatrixBCRS<float,2,2>::gpu_mult ( const TypeVector * x, 
                                                       TypeVector * y,
                                                       unsigned int vec_size) {
#ifdef USE_TEXTURE
    // store the x vector in texture
    bind_x((const float2*)x);
#endif
    const unsigned int BLOCK_SIZE = 256;
    const dim3 grid = make_large_grid(vec_size, BLOCK_SIZE);
    CNCMat2x2VecMult4Kernel<<<grid, BLOCK_SIZE>>>(
        (float4 *)(gpu_mat_), 
   	    a.size(),
	    gpu_redundant_rp_, 
	    rowptr.size(),
	    gpu_ci_, 
	    colind.size(),
	    (const float2 *)(x), 
	    (float2 *)(y), 
	    vec_size ) ;

#ifdef USE_TEXTURE
    // store the x vector in texture
    unbind_x((const float2*)x);
#endif
}

//---------------------------------------------------------------------------//
// GPU BCRS double 2x2 Matrix / Vector Multiply on GPU                       //
//---------------------------------------------------------------------------//

template <> template <typename TypeVector>
inline void CNCSparseMatrixBCRS<double,2,2>::gpu_mult ( const TypeVector * x, 
                                                        TypeVector * y,
                                                        unsigned int vec_size) {

#ifdef USE_TEXTURE
    // store the x vector in texture
    bind_x((const double2*)x);
#endif
    const unsigned int BLOCK_SIZE = 256 ;
    const dim3 grid = make_large_grid(vec_size, BLOCK_SIZE);
    CNCMat2x2VecMult4Kernel<<<grid, BLOCK_SIZE>>>(
		(double2 *)gpu_mat0_, 
		(double2 *)gpu_mat1_ , 
		a.size(),
		gpu_redundant_rp_, 
		rowptr.size(),
		gpu_ci_, 
		colind.size(),
		(const double2 *)x, 
		(double2 *)y, 
		vec_size ) ;

#ifdef USE_TEXTURE
    // store the x vector in texture
    unbind_x((const double2*)x);
#endif
}

//---------------------------------------------------------------------------//
// Conversion from NLSparseMatrix to CNCSparseMatrixBCRS : definition        //
//---------------------------------------------------------------------------//

template <class TypeMatrix, int BM, int BN> 
inline void convert_matrix( NLSparseMatrix * rhs, 
                            CNCSparseMatrixBCRS<TypeMatrix, BM, BN>& A ) {

    unsigned int BLOC_SIZE = CNCSparseMatrixBCRS<TypeMatrix, BM, BN>::BLOC_SIZE ;

    // Compute number of bloc rows and bloc columns
    unsigned int M = rhs->m / BM ;
    if((rhs->m % BM) != 0) {
        M++ ;
    }

    unsigned int N = rhs->n / BN ;
    if((rhs->n % BN) != 0) {
        N++ ;
    }

    A.N_ = N ;

    // Step 1: determine blocs to use

	CNCArray1d< std::set<unsigned int> > row_blocs(M) ;
    for(unsigned int i=0; i<rhs->m; i++) {
        unsigned int I = i / BM ;
        const NLRowColumn * Ri = &(rhs->row[i]) ;
		for(unsigned int jj=0 ; jj < Ri->size ; jj++) {
            unsigned int j = Ri->coeff[jj].index ;
            unsigned int J = j / BN ;
            row_blocs[I].insert(J) ;
        }
    }

    // Step 2: initialize rowptr 
    A.rowptr.allocate(M+1) ;
    A.rowptr[0] = 0 ;
    for(unsigned int I=0; I<M; I++) {
        A.rowptr[I+1] = (unsigned int)(A.rowptr[I] + row_blocs[I].size()) ;
    }
    unsigned int NNZ = A.rowptr[M] ;

    // Step 3: initialize colind
    A.colind.allocate(NNZ) ;
    unsigned int cur = 0 ;
    for(unsigned int I=0; I<M; I++) {
        for(std::set<unsigned int>::iterator it = row_blocs[I].begin(); it != row_blocs[I].end(); it++) {
            A.colind[cur++] = (*it) ;
        }
    }           
            
    // Step 4: initialize a
    A.a.allocate(NNZ * BLOC_SIZE) ;
    A.a.set_all(0.0) ;
    for(unsigned int i=0; i<rhs->m; i++) {
        unsigned int I = i / BM ;
        unsigned int di = i % BM ;
        const NLRowColumn * Ri = &(rhs->row[i]) ; 
		for(unsigned int jj=0; jj < Ri->size; jj++) {
            unsigned int j = Ri->coeff[jj].index ;
            unsigned int J = j / BN ;
            unsigned int dj = j % BN ;
		    for(unsigned int K=A.rowptr[I]; K<A.rowptr[I+1]; K++) {
                if(A.colind[K] == J) {
                    A.a[K * BLOC_SIZE + di * BN + dj] = static_cast<TypeMatrix>(Ri->coeff[jj].value) ;
			        break ;
                }
            }
        }
    }
    
    // Step 5: initialize diag
    A.diag.allocate(rhs->diag_size) ;
    for(unsigned int i=0; i<A.diag.size(); i++) {
        A.diag[i] = static_cast<TypeMatrix> (rhs->diag[i]) ;
    }
}

#endif
