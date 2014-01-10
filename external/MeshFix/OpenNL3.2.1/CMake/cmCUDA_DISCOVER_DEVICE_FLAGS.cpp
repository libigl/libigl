
#include <fstream>
#include <sstream>
#include <string>
#if WIN32
#include <windows.h>
#else
//#include <iostream>
#include <dirent.h>
#include <sys/stat.h>
#endif /*WIN32*/

#include "cmCPluginAPI.h"

#include <cuda.h>
#include <cuda_runtime_api.h>

class FlagBuilder{
public:
    void add_flag(const std::string& flag){
        oss << "-D" << flag << ";";
    }

    template <typename T> void add_def(const std::string& def, const T value){
        oss << "-D" << def << "=" << value << ";";
    }

    std::string str(){
        return oss.str();
    }

private:
    std::ostringstream oss;
};

FlagBuilder fb;

// Check if there is a device supporting CUDA
void GetCUDADeviceFlags()
{
    int deviceCount;
    bool archi_13 = false; 
    cudaGetDeviceCount(&deviceCount);

    // This function call returns 0 if there are no CUDA capable devices.
    if (deviceCount == 0)
        printf("There is no device supporting CUDA\n");
    int dev = 0;
    for (dev = 0; dev < 1/*deviceCount/*HACK : we want the first one only for now*/; ++dev) {
    printf("\nFirst device is selected by default.\n\n");
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);

        if (dev == 0) {
			// This function call returns 9999 for both major & minor fields, if no CUDA capable devices are present
            if (deviceProp.major == 9999 && deviceProp.minor == 9999)

            printf("There is no device supporting CUDA.\n");
            else if (deviceCount == 1)
                printf("There is 1 device supporting CUDA\n");
            else
                printf("There are %d devices supporting CUDA\n", deviceCount);
           if ((deviceProp.major * 10 + deviceProp.minor) >= 13 ) {
               // OK on peut faire du double et rajouter des optimizations
	       archi_13 = true; 
	   } 
        #if CUDART_VERSION >= 2000
        if (archi_13){
            fb.add_def("MAX_THREADS", deviceProp.multiProcessorCount*1024);
	    } else {
            fb.add_def("MAX_THREADS", deviceProp.multiProcessorCount*768);
        }
        #endif
        fb.add_def("WARP_SIZE", deviceProp.warpSize);
	fb.add_flag("USE_CRS_SHARED");
	fb.add_flag("USE_TEXTURE");
	}
        printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);
    #if CUDART_VERSION >= 2020
		int driverVersion = 0, runtimeVersion = 0;
		cudaDriverGetVersion(&driverVersion);
		printf("  CUDA Driver Version:                           %d.%d\n", driverVersion/1000, driverVersion%100);
		cudaRuntimeGetVersion(&runtimeVersion);
		printf("  CUDA Runtime Version:                          %d.%d\n", runtimeVersion/1000, runtimeVersion%100);
    #endif

        printf("  CUDA Capability Major revision number:         %d\n", deviceProp.major);
        printf("  CUDA Capability Minor revision number:         %d\n", deviceProp.minor);
        
		printf("  Total amount of global memory:                 %u bytes\n", deviceProp.totalGlobalMem);
    #if CUDART_VERSION >= 2000
        printf("  Number of multiprocessors:                     %d\n", deviceProp.multiProcessorCount);
        printf("  Number of cores:                               %d\n", 8 * deviceProp.multiProcessorCount);
    #endif
        printf("  Total amount of constant memory:               %u bytes\n", deviceProp.totalConstMem); 
        printf("  Total amount of shared memory per block:       %u bytes\n", deviceProp.sharedMemPerBlock);
        printf("  Total number of registers available per block: %d\n", deviceProp.regsPerBlock);
        printf("  Warp size:                                     %d\n", deviceProp.warpSize);
        printf("  Warp size:                                     %d\n", deviceProp.warpSize);
        printf("  Maximum number of threads per block:           %d\n", deviceProp.maxThreadsPerBlock);
        printf("  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
               deviceProp.maxThreadsDim[0],
               deviceProp.maxThreadsDim[1],
               deviceProp.maxThreadsDim[2]);
        printf("  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
               deviceProp.maxGridSize[0],
               deviceProp.maxGridSize[1],
               deviceProp.maxGridSize[2]);
        printf("  Maximum memory pitch:                          %u bytes\n", deviceProp.memPitch);
        printf("  Texture alignment:                             %u bytes\n", deviceProp.textureAlignment);
        printf("  Clock rate:                                    %.2f GHz\n", deviceProp.clockRate * 1e-6f);
    #if CUDART_VERSION >= 2000
        printf("  Concurrent copy and execution:                 %s\n", deviceProp.deviceOverlap ? "Yes" : "No");
    #endif
    #if CUDART_VERSION >= 2020
        printf("  Run time limit on kernels:                     %s\n", deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No");
        printf("  Integrated:                                    %s\n", deviceProp.integrated ? "Yes" : "No");
        printf("  Support host page-locked memory mapping:       %s\n", deviceProp.canMapHostMemory ? "Yes" : "No");
        printf("  Compute mode:                                  %s\n", deviceProp.computeMode == cudaComputeModeDefault ?
			                                                            "Default (multiple host threads can use this device simultaneously)" :
		                                                                deviceProp.computeMode == cudaComputeModeExclusive ?
																		"Exclusive (only one host thread at a time can use this device)" :
		                                                                deviceProp.computeMode == cudaComputeModeProhibited ?
																		"Prohibited (no host thread can use this device)" :
																		"Unknown");
    #endif
    }
}


#ifdef  __cplusplus
extern "C" {
#endif

static int InitialPass(void *inf, void *mf, int argc, char *argv[])
{
	cmLoadedCommandInfo *info = (cmLoadedCommandInfo *)(inf);
	const char *pwd = info->CAPI->GetCurrentDirectory(mf);
	info->CAPI->DisplaySatus(mf, pwd);
	
	GetCUDADeviceFlags();
	info->CAPI->DisplaySatus(mf, fb.str().c_str());
	info->CAPI->AddDefinition(mf, "CUDA_DEVICE_FLAGS", fb.str().c_str());
	return 1;
}


void CM_PLUGIN_EXPORT
CUDA_DISCOVER_DEVICE_FLAGSInit(cmLoadedCommandInfo *info)
{
	info->InitialPass = InitialPass;
	info->m_Inherited = 0;
	info->Name = "CUDA_DISCOVER_DEVICE_FLAGS";
}

#ifdef  __cplusplus
}
#endif
