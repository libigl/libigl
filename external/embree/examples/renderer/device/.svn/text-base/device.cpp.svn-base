// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include "device_singleray/default.h"
#include "device.h"
#include "sys/network.h"
#include "sys/library.h"

namespace embree
{
  typedef Device* (*create_device_func)(const char* parms, size_t numThreads, size_t verbose);

  Device* rtCreateDeviceHelper(const char* file, const char* parms, size_t numThreads, size_t verbose)
  {
    lib_t lib = openLibrary(file);
    if (lib == NULL) throw std::runtime_error("failed loading library \""+std::string(file)+"\"");
    create_device_func f = (create_device_func) getSymbol(lib, "create");
    if (f == NULL) throw std::runtime_error("invalid device library");
    Device* dev = f(parms,numThreads,verbose);
    if (dev == NULL) throw std::runtime_error("device creation failed");
    return dev;
  }

  Device* Device::rtCreateDevice(const char* type, size_t numThreads, size_t verbose)
  {
    if      (!strcmp(type,"default"      )) return rtCreateDeviceHelper("device_singleray","",numThreads,verbose);
    else if (!strcmp(type,"singleray"    )) return rtCreateDeviceHelper("device_singleray","",numThreads,verbose);
    else if (!strcmp(type,"singleray_knc")) return rtCreateDeviceHelper("device_coi","device_singleray_knc",numThreads,verbose);
    else if (!strcmp(type,"ivl_sse"      )) return rtCreateDeviceHelper("device_ivl_sse","",numThreads,verbose);
    else if (!strcmp(type,"ivl_avx"      )) return rtCreateDeviceHelper("device_ivl_avx","",numThreads,verbose);
    else if (!strcmp(type,"ivl_knc"      )) return rtCreateDeviceHelper("device_ivl_knc","",numThreads,verbose);
    else if (!strcmp(type,"ispc"         )) return rtCreateDeviceHelper("device_ispc","",numThreads,verbose);
    else if (!strcmp(type,"ispc_sse"     )) return rtCreateDeviceHelper("device_ispc_sse","",numThreads,verbose);
    else if (!strcmp(type,"ispc_avx"     )) return rtCreateDeviceHelper("device_ispc_avx","",numThreads,verbose);
    else if (!strcmp(type,"ispc_knc"     )) return rtCreateDeviceHelper("device_coi","device_ispc_knc",numThreads,verbose);
    else if (strstr(type,"network ")==type)  return rtCreateDeviceHelper("device_network",type+8,numThreads,verbose);
    else throw std::runtime_error("unknown device: "+std::string(type));
  }
}

