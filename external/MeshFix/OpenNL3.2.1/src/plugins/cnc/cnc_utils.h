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
 
#ifndef CNC_UTILS_H
#define CNC_UTILS_H

#include "cnc_cublas_utils.h"
#include "cnc_arrays.h"


//---------------------------------------------------------------------------//


#define CNC_ASSERT(isOK,message)                        \
  if ( !(isOK) ) {                                      \
    (void)printf("ERROR!! Assert '%s' failed\n%s\n",    \
				 #isOK, message);                       \
	return false ;                                      \
  }

#define CNC_SHORT_ASSERT(isOK,message)                  \
  if ( !(isOK) ) {                                      \
    (void)printf("%s\n", message);                      \
	return false ;                                      \
  }

#define CNC_SIMPLE_ASSERT(isOK,message)                 \
  if ( !(isOK) ) {                                      \
    (void)printf("ERROR!! Assert '%s' failed\n%s\n",    \
				 #isOK, message);                       \
  }

//---------------------------------------------------------------------------//
// Memory allocation functions host / device                                 //
//---------------------------------------------------------------------------//
template <typename TypeVector>
inline void new_device_vector(TypeVector ** device_vec, unsigned int size_vec){
    cublasAlloc( size_vec , sizeof(TypeVector), (void **)device_vec ) ;
} 

template <typename TypeVector>
inline void new_host_vector(TypeVector ** host_vec, unsigned int size_vec){
    *host_vec = NL_NEW_ARRAY(TypeVector, size_vec);
    NL_CLEAR_ARRAY(TypeVector , *host_vec , size_vec);  
} 

//---------------------------------------------------------------------------//
// Memory deallocation functions host / device                               //
//---------------------------------------------------------------------------//

template <typename TypeVector>
inline void delete_device_vector(TypeVector ** device_vec){
    cublasFree( *device_vec ) ;
    *device_vec = NULL;
}

template <typename TypeVector>
inline void delete_host_vector(TypeVector ** host_vec){
    NL_DELETE_ARRAY( *host_vec);
} 

//---------------------------------------------------------------------------//
// Memory copying functions host to device                                   //
//---------------------------------------------------------------------------//

template <typename TypeVector>
inline void copy_host_to_device_vector(
        const TypeVector * host_vec,
        TypeVector * device_vec, 
        unsigned int size_vec){
    cublasSetVector( size_vec, sizeof(TypeVector), host_vec, 1, (void*)device_vec, 1 ) ;
} 

template <typename TypeVector>
inline void copy_host_to_device_vector(
        const CNCArray1d<TypeVector> & host_array,
        TypeVector * device_vec){
    cublasSetVector( host_array.size(), sizeof(TypeVector), host_array.data(), 1, (void*)device_vec, 1 ) ;
} 

//---------------------------------------------------------------------------//
// Memory copying functions device to host                                   //
//---------------------------------------------------------------------------//

template <typename TypeVector>
inline void copy_device_to_host_vector(
        const TypeVector * device_vec,
        TypeVector * host_vec,
        unsigned int size_vec){
    cublasGetVector( size_vec, sizeof(TypeVector), (void*)device_vec, 1, host_vec, 1 );
} 

template <typename TypeVector>
inline void copy_device_to_host_vector(
        const TypeVector * device_vec,
        CNCArray1d<TypeVector> & host_array){
    cublasGetVector( host_array.size(), sizeof(TypeVector), (void*)device_vec, 1, host_array.data(), 1 );
}

//---------------------------------------------------------------------------//
// Wrappers for cublas : wrapXfunc -> cublas + {S|D} + {axpy | scal | dot }
//---------------------------------------------------------------------------//

// cublasSaxpy
inline void wrapXaxpy ( int N, float alpha, const float * gpu_x, int inc_x, float * gpu_y, int inc_y){
    cublasSaxpy( N, alpha, gpu_x, inc_x, gpu_y, inc_y);
}

// cublasDaxpy
inline void wrapXaxpy ( int N, double alpha, const double * gpu_x, int inc_x, double * gpu_y, int inc_y){
    cublasDaxpy( N, alpha, gpu_x, inc_x, gpu_y, inc_y);
}

// cublasScal
inline void wrapXscal ( int N, float alpha, float * gpu_x, int inc_x){
    cublasSscal( N, alpha, gpu_x, inc_x);
}

// cublasDscal
inline void wrapXscal ( int N, double alpha, double * gpu_x, int inc_x){
    cublasDscal( N, alpha, gpu_x, inc_x);
}

// cublasSdot
inline float wrapXdot ( int N, const float * x, int inc_x, const float *y, int inc_y){
    return  cublasSdot( N, x, inc_x ,y, inc_y);
}

// cublasDdot
inline double wrapXdot ( int N, const double * x, int inc_x, const double *y, int inc_y){
    return  cublasDdot( N, x, inc_x ,y, inc_y);
}

//---------------------------------//
// Wrapper for the sqrt function  //
//-------------------------------- // 
inline float sqrt_gen ( float x){
    return  sqrtf(x);
}

inline double sqrt_gen ( double x){
    return  sqrt(x);
}


#endif
