// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#include "platform.h"
#include "intrinsics.h"

////////////////////////////////////////////////////////////////////////////////
/// Windows Platform
////////////////////////////////////////////////////////////////////////////////

#ifdef _WIN32

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

namespace embree
{
  void* os_malloc(size_t bytes) {
    char* ptr = (char*) VirtualAlloc(NULL,bytes,MEM_COMMIT|MEM_RESERVE,PAGE_READWRITE);
    if (ptr == NULL) throw std::runtime_error("memory allocation failed");
    for (size_t i=0; i<bytes; i+=4096) ptr[i] = 0; // touch pages
    return ptr;
  }

  void os_free(void* ptr, size_t bytes) {
    VirtualFree(ptr,bytes,MEM_RELEASE);
  }

  double getSeconds() {
    LARGE_INTEGER freq, val;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&val);
    return (double)val.QuadPart / (double)freq.QuadPart;
  }
}
#endif

////////////////////////////////////////////////////////////////////////////////
/// Unix Platform
////////////////////////////////////////////////////////////////////////////////

#if defined(__UNIX__)

#include <sys/time.h>
#include <sys/mman.h>

namespace embree
{
  /*! page allocation */
  void* os_malloc(size_t bytes) {
    char* ptr = (char*) mmap(0, bytes, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANON, -1, 0);
    if (ptr == NULL) throw std::runtime_error("memory allocation failed");
    for (size_t i=0; i<bytes; i+=4096) ptr[i] = 0; // touch pages                                                                             
    return ptr;
  }

  void os_free(void* ptr, size_t bytes) {
    munmap(ptr,bytes);
  }

  double getSeconds() {
    struct timeval tp; gettimeofday(&tp,NULL);
    return double(tp.tv_sec) + double(tp.tv_usec)/1E6;
  }
}

#endif

////////////////////////////////////////////////////////////////////////////////
/// All Platforms
////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <string.h>

namespace embree
{
  void* alignedMalloc(size_t size, size_t align)
  {
    char* base = (char*)malloc(size+align+sizeof(int));
    if (base == NULL) throw std::runtime_error("memory allocation failed");
    char* unaligned = base + sizeof(int);
    char*   aligned = unaligned + align - ((size_t)unaligned & (align-1));
    ((int*)aligned)[-1] = (int)((size_t)aligned-(size_t)base);
    return aligned;
  }
  
  void alignedFree(const void* ptr) {
    int ofs = ((int*)ptr)[-1];
    free((char*)ptr-ofs);
  }
}
