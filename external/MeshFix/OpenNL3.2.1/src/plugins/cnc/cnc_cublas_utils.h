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
#ifndef CNC_GPU_UTILS_H
#define CNC_GPU_UTILS_H

#include <cublas.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <string.h>


//---------------------------------------------------------------------------//
// thread block size (number of threads per block = BLOCK_SIZE*BLOCK_SIZE)	 //
//---------------------------------------------------------------------------//

#define THREAD_BLOCK_SIZE 16


static void CNCgetError ( cublasStatus st ) {
	switch (st) {
		case CUBLAS_STATUS_SUCCESS:	break ;  //printf ( "cublas : no error\n" ) ; break ;
		case CUBLAS_STATUS_NOT_INITIALIZED:	printf ("CUBLAS_STATUS_NOT_INITIALIZED\n") ;break ;
		case CUBLAS_STATUS_ALLOC_FAILED:	printf ("CUBLAS_STATUS_ALLOC_FAILED\n") ;	break ;
		case CUBLAS_STATUS_INVALID_VALUE:	printf ("CUBLAS_STATUS_INVALID_VALUE\n") ;	break ;
		case CUBLAS_STATUS_MAPPING_ERROR:	printf ("CUBLAS_STATUS_MAPPING_ERROR\n") ;	break ;
		case CUBLAS_STATUS_EXECUTION_FAILED:printf ("CUBLAS_STATUS_EXECUTION_FAILED\n");break ;
		case CUBLAS_STATUS_INTERNAL_ERROR:	printf ("CUBLAS_STATUS_INTERNAL_ERROR\n") ;	break ;
	    default: printf ("unkown error message\n"); break ;
	}
}


//static void CNCgetError () {
//	cublasStatus st = cublasGetError() ;
//	CNCgetError ( st ) ;
//}


// Check if there is a device supporting CUDA
static bool CNCCheckDevice()
{
#if __DEVICE_EMULATION__
    return true;
#else
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0){
        return false;
    }
    int dev;
    for (dev = 0; dev < deviceCount; ++dev) {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);
        if (strncmp(deviceProp.name, "Device Emulation", 16))
            break;
    }
    if (dev == deviceCount) {
        return false;
    } else {
        cudaSetDevice(dev);
        return true;
    }
#endif
}

// Check if the device that would be found by CNCCheckDevice supports doubles
static bool CNCCheckDeviceDoubleSupport()
{
#if __DEVICE_EMULATION__
    return true;
#else
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0){
        return false;
    }
    int dev;
    for (dev = 0; dev < deviceCount; ++dev) {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);
        if (strncmp(deviceProp.name, "Device Emulation", 16)){
            if(deviceProp.major > 1 || deviceProp.minor >= 3){
                return true;
            }
            return false;
        }
    }
    return false;
#endif
}

// dummy kernel used to setup L1 cache on Fermi cards
#if CUDA_VERSION >= 3000
extern "C" __global__ static void dummy_kernel() { }
#endif

// Configure a CUDA device to have the best performance
static bool CNCConfigureDevice()
{
#if __DEVICE_EMULATION__
    return true;
#else
    int deviceCount, dev, devCurrent;
    cudaDeviceProp deviceProp;

    // Save current device
    cudaGetDevice(&devCurrent);

    // Loop through all the available devices
    cudaGetDeviceCount(&deviceCount);
    for (dev=0; dev < deviceCount; ++dev) {
        cudaGetDeviceProperties(&deviceProp, dev);

        // Use 48k of L1 cache on Fermi cards. We only set this on a dummy
        // kernel, but this setting is persistent and will thus affect all our
        // kernels. Launching this kernel is required for this setting to be
        // taken into account.
#if CUDA_VERSION >= 3000
        if (deviceProp.major >= 2) {
            cudaSetDevice(dev);
            if (cudaFuncSetCacheConfig(dummy_kernel, cudaFuncCachePreferL1) != cudaSuccess)
                return false;
            dummy_kernel<<< 8192, 256 >>>();
            if (cudaThreadSynchronize() != cudaSuccess)
                return false;
        }
    }
#endif

    // Restore current device
    cudaSetDevice(devCurrent);

    return true;
#endif
}



#endif //CNC_UTILS_H
