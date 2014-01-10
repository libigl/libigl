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



#ifndef CNC_TEXTURE_H
#define CNC_TEXTURE_H

// taken from Nathan Bell's spmv 

#ifdef __CUDACC__
#include <cuda.h>
#  define CUDA_SAFE_CALL_NO_SYNC( call) do {                                 \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } } while (0)

#  define CUDA_SAFE_CALL( call) do {                                         \
    CUDA_SAFE_CALL_NO_SYNC(call);                                            \
    cudaError err = cudaThreadSynchronize();                                 \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } } while (0)

#else
// drop all CUDA calls
#define CUDA_SAFE_CALL_NO_SYNC(x)
#define CUDA_SAFE_CALL(x)
#endif

texture<float,1> tex_x_float;
texture<float4,1> tex_x_float4;
texture<float2,1> tex_x_float2;
// Use int2 to pull doubles through texture cache
texture<int2,1>  tex_x_double;
// Use int4 to pull double2s through texture cache
texture<int4,1> tex_x_double2;


void bind_x(const float * x) {   
   CUDA_SAFE_CALL(cudaBindTexture(NULL, tex_x_float, x));   
}
void bind_x(const float2 * x) {   
   CUDA_SAFE_CALL(cudaBindTexture(NULL, tex_x_float2, x));   
}
void bind_x(const float4 * x) {   
   CUDA_SAFE_CALL(cudaBindTexture(NULL, tex_x_float4, x));   
}

void bind_x(const double * x) {   
   CUDA_SAFE_CALL(cudaBindTexture(NULL, tex_x_double, x));   
} 

void bind_x(const double2 * x) {   
   CUDA_SAFE_CALL(cudaBindTexture(NULL, tex_x_double2, x));   
} 

void unbind_x(const float * x) { 
   CUDA_SAFE_CALL(cudaUnbindTexture(tex_x_float)); 
}

void unbind_x(const float4 * x) { 
   CUDA_SAFE_CALL(cudaUnbindTexture(tex_x_float4)); 
}

void unbind_x(const float2 * x) { 
   CUDA_SAFE_CALL(cudaUnbindTexture(tex_x_float2)); 
}

void unbind_x(const double * x) {   
   CUDA_SAFE_CALL(cudaUnbindTexture(tex_x_double)); 
}
void unbind_x(const double2 * x) {   
   CUDA_SAFE_CALL(cudaUnbindTexture(tex_x_double2)); 
}

__inline__ __device__ float fetch_x(const int& i, const float * x) {
   return tex1Dfetch(tex_x_float, i);
}

__inline__ __device__ float4 fetch_x(const int& i, const float4 * x) {
   return tex1Dfetch(tex_x_float4, i);
}

__inline__ __device__ float2 fetch_x(const int& i, const float2 * x) {
   return tex1Dfetch(tex_x_float2, i);
}

#if !defined(CUDA_NO_SM_13_DOUBLE_INTRINSICS)
__inline__ __device__ double fetch_x(const int& i, const double * x)
{
    int2 v = tex1Dfetch(tex_x_double, i);
    return __hiloint2double(v.y, v.x);
}

__inline__ __device__ double2 fetch_x(const int& i, const double2 * x)
{
    int4 v=tex1Dfetch(tex_x_double2, i);
	return make_double2(__hiloint2double(v.y, v.x),__hiloint2double(v.w, v.z));
}
#else // !defined(CUDA_NO_SM_13_DOUBLE_INTRINSICS)

// dummy implementation just for compiling 
__inline__ __device__ double fetch_x(const int& i, const double * x)
{
    return 0.;
}

__inline__ __device__ double2 fetch_x(const int& i, const double2 * x)
{
	return make_double2(0.,0.);
}
#endif



#endif
