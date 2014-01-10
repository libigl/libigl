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


#ifndef CNC_KERNELS_H
#define CNC_KERNELS_H

#include <functional>
#include <algorithm>
#include "cnc_texture.h"


#define THREAD_BLOCK_SIZE 16
#define CUDA_BLOCK_SIZE 256
#ifndef WARP_SIZE
  #define WARP_SIZE 32 
#endif
#ifndef MAX_THREADS
  #define MAX_THREADS (4*768)
#endif
//#define USE_SHARED_CRS   
//#define USE_TEXTURE

#ifdef __DEVICE_EMULATION__
#define EMUSYNC __syncthreads()
#else   
#define EMUSYNC
#endif




//#define small_grid_thread_id(void) ((blockDim.x * blockIdx.x + threadIdx.x))
//#define large_grid_thread_id(void) ((blockDim.x * (blockIdx.x + blockIdx.y*gridDim.x) + threadIdx.x))
#define small_grid_thread_id(void) ((__umul24(blockDim.x, blockIdx.x) + threadIdx.x))
#define large_grid_thread_id(void) ((__umul24(blockDim.x,blockIdx.x + __umul24(blockIdx.y,gridDim.x)) + threadIdx.x))
#define DIVIDE_INTO(x,y) ((x + y - 1)/y)



//---------------------------------------------------------------------------//
// to make a grid a block : taken from Nathan Bell  spmv program             //
//---------------------------------------------------------------------------//

dim3 make_large_grid(const unsigned int num_threads, const unsigned int blocksize){
    const unsigned int num_blocks = DIVIDE_INTO(num_threads, blocksize);
    if (num_blocks <= 65535){
        //fits in a 1D grid
        return dim3(num_blocks);
    } else {
        //2D grid is required
        const unsigned int side = (unsigned int) ceil(sqrt((double)num_blocks));
        return dim3(side,side);
    }
}

//---------------------------------------------------------------------------//
// kernel matrix-vector product for FLOAT BCRS2x2 matrix                     //
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
        unsigned int size_vec ); 

//---------------------------------------------------------------------------//
// kernel matrix-vector product for DOUBLE BCRS2x2 matrix                    //
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
        unsigned int size_vec );


//---------------------------------------------------------------------------//
// template kernel for matrix-vector product with a CRS matrix               //       
//---------------------------------------------------------------------------//
template <typename TypeMatrix,typename TypeVector> 
__global__ void CNCMat1x1VecMultKernel (
        TypeMatrix * matrix,
        unsigned int size_matrix,
        uint2 * rowptr,
        unsigned int size_rowptr,
        unsigned int * colind,
        unsigned int size_colind,
        const TypeVector * x,
        TypeVector * b, 
        unsigned int size_vec) {

	// Thread index
	const unsigned int index = large_grid_thread_id(void);

	if ( index < size_vec ) {

		uint2 rowptr_bounds = rowptr[index] ;

		TypeVector res = 0.;

		// for each block of the block_row, mult
		for ( unsigned int i=rowptr_bounds.x; i<rowptr_bounds.y; i++ ) { 
            #ifndef USE_TEXTURE
	    	res += static_cast<TypeVector>(matrix[i]) * x[colind[i]] ;
            #else
			res += static_cast<TypeVector>(matrix[i]) * fetch_x(colind[i],x) ;
            #endif
		}
		b[index] = res ;
	}
}

//---------------------------------------------------------------------------//
// template kernel vector for matrix-product vector with a CRS matrix        //       
// taken from Nathan Bell spmv program                                       //
//---------------------------------------------------------------------------//
template <typename TypeMatrix,typename TypeVector,unsigned int BLOCK_SIZE> 
__global__ void CNCCRSNathan_Bell (
        TypeMatrix * matrix,
        unsigned int size_matrix,
        uint2 * rowptr,
        unsigned int size_rowptr,
        unsigned int * colind,
        unsigned int size_colind,
        const TypeVector * x, 
        TypeVector * b, 
        unsigned int size_vec) {

    __shared__ TypeMatrix  sdata[BLOCK_SIZE+16];
    __shared__ uint2 ptrs[BLOCK_SIZE/WARP_SIZE];
    
    const unsigned int thread_id   = BLOCK_SIZE * blockIdx.x + threadIdx.x;  // global thread index
    const unsigned int thread_lane = threadIdx.x & (WARP_SIZE-1);            // thread index within the warp
    const unsigned int warp_id     = thread_id   / WARP_SIZE;                // global warp index
    const unsigned int warp_lane   = threadIdx.x / WARP_SIZE;                // warp index within the CTA
    const unsigned int num_warps   = (BLOCK_SIZE / WARP_SIZE) * gridDim.x;   // total number of active warps

    for(unsigned int index = warp_id; index < size_vec; index += num_warps){
        // use two threads to fetch Ap[index] and Ap[index+1]
        // this is considerably faster than the more straightforward option
        if (thread_lane < 2){
	        ptrs[warp_lane] = rowptr[index];
        }
        
        const unsigned int row_start = ptrs[warp_lane].x; //same as: row_start = Ap[index];
        const unsigned int row_end   = ptrs[warp_lane].y; //same as: row_end   = Ap[index+1];

        // compute local sum
        TypeVector sum = 0;
        for(unsigned int  jj = row_start + thread_lane; jj < row_end; jj += WARP_SIZE){
            #ifndef USE_TEXTURE
            sum += static_cast<TypeVector>(matrix[jj]) * x[colind[jj]];
            #else
            sum += static_cast<TypeVector>(matrix[jj]) * fetch_x(colind[jj],x);
            #endif
        }
        
        // reduce local sums to row sum (ASSUME: warpsize 32)
	    sdata[threadIdx.x] = sum;
        sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x + 16]; EMUSYNC; 
        sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  8]; EMUSYNC; 
        sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  4]; EMUSYNC; 
        sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  2]; EMUSYNC; 
        sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  1]; EMUSYNC; 

        // first thread writes warp result
        if (thread_lane == 0){
            b[index] = sdata[threadIdx.x];
        }
    }
}


//---------------------------------------------------------------------------//
// template kernel vector for matrix-product vector with a ELL matrix        //       
// taken from Nathan Bell article                                            //
//---------------------------------------------------------------------------//
template <typename TypeMatrix,typename TypeVector> 
__global__ void CNCMatELL_Nathan_Bell (
        const unsigned int num_rows, 
        const unsigned int num_cols,
        const unsigned int num_cols_per_row, 
        const unsigned int stride,
        TypeMatrix * matrix,
        unsigned int * indices, 
        const TypeVector * x,
        TypeVector * y) {

    const unsigned int row = large_grid_thread_id(void);
    if (row < num_rows) {
       
        TypeVector accu = 0;
       
        for(unsigned int j = 0 ; j < num_cols_per_row ; ++j) {
           
            unsigned int index = indices[ stride * j + row ];  
      	    TypeVector   val   =  matrix[ stride * j + row ];
	        if ( val != 0) {
                  #ifndef USE_TEXTURE
                  accu += static_cast<TypeVector>(val) * x[index] ;
                  #else
                  accu += static_cast<TypeVector>(val) * fetch_x(index, x);
                  #endif
            }
        }
        y[row] = accu;
    }
}

//---------------------------------------------------------------------------//
// segmented reduction in shared memory for COO kernel                       //       
// taken from Nathan Bell                                                    //
//---------------------------------------------------------------------------//

template <typename VectorType>
__device__ VectorType segreduce_warp(
        const unsigned int thread_lane,
        unsigned int row,
        VectorType val,
        unsigned int * rows,
        VectorType * vals){
    rows[threadIdx.x] = row;
    vals[threadIdx.x] = val;

    if( thread_lane >=  1 && row == rows[threadIdx.x -  1] ) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  1]; } 
    if( thread_lane >=  2 && row == rows[threadIdx.x -  2] ) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  2]; }
    if( thread_lane >=  4 && row == rows[threadIdx.x -  4] ) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  4]; }
    if( thread_lane >=  8 && row == rows[threadIdx.x -  8] ) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  8]; }
    if( thread_lane >= 16 && row == rows[threadIdx.x - 16] ) { vals[threadIdx.x] = val = val + vals[threadIdx.x - 16]; }

    return val;
}

//---------------------------------------------------------------------------//
// subroutine for the COO kernel                                             //       
// taken from Nathan Bell                                                    //
//---------------------------------------------------------------------------//

template <typename VectorType>
__device__ void segreduce_block(const unsigned int * idx, VectorType * val)
{
    VectorType left = 0;
    if( threadIdx.x >=   1 && idx[threadIdx.x] == idx[threadIdx.x -   1] ) { left = val[threadIdx.x -   1]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();  
    if( threadIdx.x >=   2 && idx[threadIdx.x] == idx[threadIdx.x -   2] ) { left = val[threadIdx.x -   2]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >=   4 && idx[threadIdx.x] == idx[threadIdx.x -   4] ) { left = val[threadIdx.x -   4]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >=   8 && idx[threadIdx.x] == idx[threadIdx.x -   8] ) { left = val[threadIdx.x -   8]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >=  16 && idx[threadIdx.x] == idx[threadIdx.x -  16] ) { left = val[threadIdx.x -  16]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >=  32 && idx[threadIdx.x] == idx[threadIdx.x -  32] ) { left = val[threadIdx.x -  32]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();  
    if( threadIdx.x >=  64 && idx[threadIdx.x] == idx[threadIdx.x -  64] ) { left = val[threadIdx.x -  64]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >= 128 && idx[threadIdx.x] == idx[threadIdx.x - 128] ) { left = val[threadIdx.x - 128]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >= 256 && idx[threadIdx.x] == idx[threadIdx.x - 256] ) { left = val[threadIdx.x - 256]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
}

//---------------------------------------------------------------------------//
// template kernel vector for matrix-product vector with a COO matrix        //       
// taken from Nathan Bell spmv                                               //
//---------------------------------------------------------------------------//
template <typename TypeMatrix, typename TypeVector, unsigned int BLOCK_SIZE>
__global__ void spmv_coo_flat_kernel(
        const unsigned int num_nonzeros,
        const unsigned int interval_size,
        const unsigned int * I, 
        const unsigned int * J, 
        const TypeMatrix * V, 
        const TypeVector * x, 
        TypeVector * y,
        unsigned int * temp_rows,
        TypeVector * temp_vals){
        
    __shared__ unsigned int rows[BLOCK_SIZE];
    __shared__ TypeVector vals[BLOCK_SIZE];

    const unsigned int thread_id   = BLOCK_SIZE * blockIdx.x + threadIdx.x;                // global thread index
    const unsigned int thread_lane = threadIdx.x & (WARP_SIZE-1);                          // thread index within the warp
    const unsigned int warp_id     = thread_id   / WARP_SIZE;                              // global warp index

    const unsigned int interval_begin = warp_id * interval_size;                           // warp's offset into I,J,V
    const unsigned int interval_end   = min(interval_begin + interval_size, num_nonzeros); // end of warps's work

    if(interval_begin >= interval_end){                                                     // warp has no work to do 
        return;
    }

    if (thread_lane == 31){
        // initialize the carry in values
        rows[threadIdx.x] = I[interval_begin]; 
        vals[threadIdx.x] = 0;
    }
  
    for(unsigned int n = interval_begin + thread_lane; n < interval_end; n += WARP_SIZE){
        unsigned int row = I[n];                                         // row index (i)
        #ifndef USE_TEXTURE
	    TypeVector val = static_cast<TypeVector>(V[n]) * x[J[n]];            // A(i,j) * x(j)
        #else 
	    TypeVector val = static_cast<TypeVector>(V[n]) * fetch_x(J[n], x);            // A(i,j) * x(j)
        #endif 
        if (thread_lane == 0){
            if(row == rows[threadIdx.x + 31]){
                val += vals[threadIdx.x + 31];                        // row continues
            } else {
                y[rows[threadIdx.x + 31]] += vals[threadIdx.x + 31];  // row terminated
            }
        }
        
        val = segreduce_warp(thread_lane, row, val, rows, vals);      // segmented reduction in shared memory

        if(thread_lane < 31 && row != rows[threadIdx.x + 1]){
            y[row] += val;                                            // row terminated
        }
    }

    if(thread_lane == 31){
        // write the carry out values
        temp_rows[warp_id] = rows[threadIdx.x];
        temp_vals[warp_id] = vals[threadIdx.x];
    }
}

// The second level of the segmented reduction operation

template <typename TypeVector, unsigned int BLOCK_SIZE>
__global__ void spmv_coo_reduce_update_kernel(
        const unsigned int num_warps,
        const unsigned int * temp_rows,
        const TypeVector * temp_vals,
        TypeVector * y){
    __shared__ unsigned int rows[BLOCK_SIZE + 1];    
    __shared__ TypeVector vals[BLOCK_SIZE + 1];    

    const unsigned int end = num_warps - (num_warps & (BLOCK_SIZE - 1));

    if (threadIdx.x == 0){
        rows[BLOCK_SIZE] = (unsigned int) -1;
        vals[BLOCK_SIZE] = (TypeVector)  0;
    }
    
    __syncthreads();

    unsigned int i = threadIdx.x;

    while (i < end){
        // do full blocks
        rows[threadIdx.x] = temp_rows[i];
        vals[threadIdx.x] = temp_vals[i];

        __syncthreads();

        segreduce_block(rows, vals);

        if (rows[threadIdx.x] != rows[threadIdx.x + 1]){
            y[rows[threadIdx.x]] += vals[threadIdx.x];
        }

        __syncthreads();

        i += BLOCK_SIZE; 
    }

    if (end < num_warps){
        if (i < num_warps){
            rows[threadIdx.x] = temp_rows[i];
            vals[threadIdx.x] = temp_vals[i];
        } else {
            rows[threadIdx.x] = (unsigned int) -1;
            vals[threadIdx.x] = (TypeVector)  0;
        }

        __syncthreads();
   
        segreduce_block(rows, vals);

        if (i < num_warps){
            if (rows[threadIdx.x] != rows[threadIdx.x + 1]){
                y[rows[threadIdx.x]] += vals[threadIdx.x];
            }
        }
    }
}


template <typename TypeMatrix,typename TypeVector> 
__global__ void spmv_coo_serial_kernel (
        const unsigned int num_nonzeros, 
        const unsigned int * I, 
        const unsigned int * J, 
        const TypeMatrix   * V, 
        const TypeVector   * x, 
        TypeVector         * y){
    for(unsigned int n = 0; n < num_nonzeros; n++){
        y[I[n]] += static_cast<TypeVector>(V[n]) * x[J[n]];
    }
}

//---------------------------------------------------------------------------//
// template kernel for Hadamard product                                      //
//---------------------------------------------------------------------------//
template <typename TypeVector> 
__global__ void CNCVecVecMultKernel (
        unsigned int size,
        const TypeVector * x,
        const TypeVector * y,
        TypeVector * r ) {

	// Thread index
	const unsigned int index = large_grid_thread_id(void) ;
		
	if ( index < size ){
		r[index] = x[index]*y[index] ;
    }
}

//---------------------------------------------------------------------------//
// template kernel for for memory device filling  with a constant value      //
//---------------------------------------------------------------------------//
template <typename TypeVector> 
__global__ void CNCfill (
        TypeVector * dest,
        unsigned int vec_size,
        const TypeVector val ) {

	// Thread index
	const unsigned int index = large_grid_thread_id(void);  // global thread index 
	if ( index < vec_size ){
        dest[index] = val;
    }
}

#endif

