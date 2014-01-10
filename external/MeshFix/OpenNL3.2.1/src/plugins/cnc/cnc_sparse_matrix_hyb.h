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
#ifndef CNC_SPARSE_MATRIX_HYB_H
#define CNC_SPARSE_MATRIX_HYB_H

#include "cnc_utils.h"
#include "cnc_cublas_utils.h"
#include "cnc_arrays.h"
#include "cnc_kernels.h"
#include "cnc_texture.h"
#include "cnc_sparse_matrix_ell.h"
#include "cnc_sparse_matrix_coo.h"
#include <ostream>
#include <iostream>
#include <algorithm>
#include <iterator>

extern "C" {
#include <NL/nl_matrix.h>
}


class CNCSparseMatrixPatternHYB {
public:  
    typedef unsigned int index_type;
    typedef CNCArray1d<index_type> IndexArray ;

    CNCSparseMatrixPatternHYB() {
        num_rows          = 0 ;
        num_cols          = 0 ;
        num_nonzeros      = 0 ;
    }

    void clear() {
        num_rows          = 0 ;
        num_cols          = 0 ;
        num_nonzeros      = 0 ;
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
    
};

template<typename TypeMatrix> 
class CNCSparseMatrixHYB : public CNCSparseMatrixPatternHYB {
public:

    typedef TypeMatrix   type_value;
    typedef CNCArray1d<TypeMatrix> CoeffArray ;

    inline CNCSparseMatrixHYB()  { }
       
    template <typename TypeVector>
    inline void mult(const TypeVector* x, 
                     TypeVector* y,
                     const unsigned int vec_size) const ;
        
    template <typename TypeVector>
    inline void mult(const CNCArray1d<TypeVector>& x, 
                     CNCArray1d<TypeVector>& y) const {
        mult(x.data(), y.data(), x.size()) ;
    }

    inline unsigned int  data_size() const {
        return ell.a.size() + coo.V.size() ;
    }    
        
    inline TypeMatrix * diagonal() {
       return ell.diag.data() ;  
    }
        
    inline const TypeMatrix * diagonal() const {
       return ell.diag.data() ;  
    }
    
    inline void gpu_allocate_and_upload () ;
    inline void gpu_deallocate () ;
        
    
    template <typename TypeVector>
    inline void gpu_mult ( const TypeVector * x, 
                           TypeVector* y, 
                           unsigned int vec_size) ;


    //inline void print(std::ostream& out) ;
    inline void print (); 

    CNCSparseMatrixELL<TypeMatrix> ell;
    CNCSparseMatrixCOO<TypeMatrix> coo;

} ;

//---------------------------------------------------------------------------//
// printf function for HYB sparse matrix                                     //
//---------------------------------------------------------------------------//


template <typename TypeMatrix> 
inline void CNCSparseMatrixHYB<TypeMatrix>::print () {
     std::cout << "ELL Part  " << std::endl;
     ell.print();
     std::cout << "COO Part  " << std::endl;
     coo.print();
}

//---------------------------------------------------------------------------//
// Matrix-Vector product for HYB sparse matrix on CPU                        //
//---------------------------------------------------------------------------//


template <typename TypeMatrix> 
template <typename TypeVector>
inline void CNCSparseMatrixHYB<TypeMatrix>::mult(
        const TypeVector* x, 
        TypeVector* y,
        const unsigned int vec_size) const {

     ell.mult( x, y, vec_size );
     TypeVector * tmp ;
     new_host_vector( &tmp, vec_size ); 
     coo.mult( x, tmp, vec_size );
     std::transform(y , y + vec_size, tmp , y , std::plus<TypeVector>());   
     delete_host_vector( &tmp );
}

//---------------------------------------------------------------------------//
// Allocate and upload a HYB sparse matrix on GPU memory                     //
//---------------------------------------------------------------------------//

template <typename TypeMatrix>
inline void CNCSparseMatrixHYB<TypeMatrix>::gpu_allocate_and_upload () {
      ell.gpu_allocate_and_upload();
      coo.gpu_allocate_and_upload();
}

//---------------------------------------------------------------------------//
// Free graphics memory used to store a HYB sparse matrix                     //
//---------------------------------------------------------------------------//

template <typename TypeMatrix> 
inline void CNCSparseMatrixHYB<TypeMatrix>::gpu_deallocate () {
      ell.gpu_deallocate();
      coo.gpu_deallocate();
}

//---------------------------------------------------------------------------//
// GPU HYB Matrix / Vector Multiply                                      // 
//---------------------------------------------------------------------------// 
template <typename TypeMatrix>
template <typename TypeVector> 
inline void CNCSparseMatrixHYB<TypeMatrix>::gpu_mult (
        const TypeVector * x, 
        TypeVector * y,
        unsigned int vec_size) {

    ell.gpu_mult( x, y, vec_size );
     
    // slow, should be done only once in gpu_allocate_and_upload
    // but TypeVector is not a template parameter of the class
    TypeVector * coo_res;
    new_device_vector( &coo_res, vec_size);
    
    coo.gpu_mult( x, coo_res, vec_size );
     
    wrapXaxpy( vec_size , 1. , coo_res , 1 , y , 1 );    
    delete_device_vector( &coo_res ); 
}

//---------------------------------------------------------------------------//
//  Heuristic to choose the number of columns per row for the ELL part       //    
//  Taken From Nathan Bell programs                                          //    
//---------------------------------------------------------------------------//

unsigned int compute_hyb_cols_per_row(const NLSparseMatrix *rhs, 
                                      float relative_speed = 3.0, 
                                      unsigned int breakeven_threshold = 4096){
    
    // compute maximum row length
    unsigned int max_cols_per_row = 0;
    for(unsigned int i = 0; i < rhs->m; i++){
        max_cols_per_row = std::max(max_cols_per_row, rhs->row[i].size ); 
    }

    // compute distribution of nnz per row
    unsigned int * histogram ;
    new_host_vector( &histogram, max_cols_per_row + 1);
    std::fill(histogram, histogram + max_cols_per_row + 1, 0);
    for(unsigned int i = 0; i < rhs->m; i++){
        histogram[rhs->row[i].size]++;
    }

    // compute optimal ELL col
    unsigned int num_cols_per_row = max_cols_per_row;
    for(unsigned int i = 0, rows = rhs->m; i < max_cols_per_row; i++){
        rows -= histogram[i];  //number of rows of length > i
        if(relative_speed * rows < rhs->m || rows < breakeven_threshold){
            num_cols_per_row = i;
            break;
        }
    }

    delete_host_vector( &histogram );

    return num_cols_per_row;
}


//---------------------------------------------------------------------------//
//  Conversion NLSparseMatrix to CNCSparseMatrixHYB<TypeMatrix>              //    
//---------------------------------------------------------------------------//

template <typename TypeMatrix> 
inline void convert_matrix( NLSparseMatrix * rhs, 
                            CNCSparseMatrixHYB<TypeMatrix>& A, 
                            const unsigned int alignment = 32 ) {

    CNC_SIMPLE_ASSERT ( rhs->n == rhs->m, "In CNC, matrices must be square!" ) ;
    
    // references on the target matrix A
    CNCSparseMatrixELL<TypeMatrix> & ell = A.ell;
    CNCSparseMatrixCOO<TypeMatrix> & coo = A.coo;
    
    // assign num_cols
    A.num_cols = rhs->n;
    ell.num_cols = rhs->n;
    coo.num_cols = rhs->n;
    
    // assign num_rows
    A.num_rows = rhs->m ;
    ell.num_rows = rhs->m ;
    coo.num_rows = rhs->m ;

    // assign nonzeros
    unsigned int num_nonzeros = nlSparseMatrixNNZ(rhs);
    A.num_nonzeros = num_nonzeros;
        
    // compute the nnz number in the ELL and COO portion
    ell.num_cols_per_row = compute_hyb_cols_per_row(rhs); 
        
    ell.num_nonzeros = 0;
    for(unsigned int i=0; i<rhs->m; i++) {
        ell.num_nonzeros += std::min(ell.num_cols_per_row,rhs->row[i].size); 
    }
      
    // we store the rest part of nnz in COO format
    coo.num_nonzeros = num_nonzeros - ell.num_nonzeros;
        
    if (coo.num_nonzeros > 0) {
        coo.I.allocate( coo.num_nonzeros );
        coo.J.allocate( coo.num_nonzeros );
        coo.V.allocate( coo.num_nonzeros );
    } else {
        // set to zeros the triplet I,J,V 
        coo.I.clear(); 
        coo.J.clear();
        coo.V.clear();
    }      

    // compute the stride 
    ell.stride = alignment * ((ell.num_rows + alignment - 1)/ alignment);
        
    unsigned int common_size = ell.stride * ell.num_cols_per_row ;
    
    // allocate the ell part
    ell.indices.allocate( common_size ) ;
    ell.a.allocate( common_size ) ;
    ell.diag.allocate( A.num_rows ) ;
       
    // filling to zero the ell part
    std::fill( ell.indices.data(), ell.indices.data() + common_size, 0);
    std::fill(       ell.a.data(),       ell.a.data() + common_size, static_cast<TypeMatrix>(0));
    std::fill(    ell.diag.data(),    ell.diag.data() + ell.num_rows, static_cast<TypeMatrix>(0));
    
    
    for(unsigned int i=0, coo_nnz = 0; i<rhs->m; i++) {
        
        NLRowColumn * R = &(rhs->row[i]) ;
        unsigned int jj = 0;
        while( jj<R->size && jj < ell.num_cols_per_row ) {
                 
            // in of (i,jj) in A
            unsigned int index = ell.stride* jj + i;
         
            ell.a[index]       = static_cast<TypeMatrix>(R->coeff[jj].value) ;
            ell.indices[index] =                         R->coeff[jj].index  ;
            // if diagonal element we store it in diag
            if (i == R->coeff[jj].index){
                ell.diag[i] = static_cast<TypeMatrix>(R->coeff[jj].value);
            }
            jj++;
        }
        while( jj < R->size) {
            coo.I[coo_nnz] = i;
            coo.J[coo_nnz] = R->coeff[jj].index;
            coo.V[coo_nnz] = static_cast<TypeMatrix>(R->coeff[jj].value) ;
            jj++;
            coo_nnz++;
        }
    }
    std::cout << "ELL num cols per rows = "  << ell.num_cols_per_row << std::endl;
    std::cout << "ELL nnz          = "       << ell.num_nonzeros << std::endl;
    std::cout << "COO nnz per rows = "       << coo.num_nonzeros << std::endl;
    //A.print();
}

#endif
