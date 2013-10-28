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
  void* os_malloc(size_t bytes) 
  {
    char* ptr = (char*) VirtualAlloc(NULL,bytes,MEM_COMMIT|MEM_RESERVE,PAGE_READWRITE);
    if (ptr == NULL) throw std::runtime_error("memory allocation failed");
    return ptr;
  }

  void* os_reserve(size_t bytes)
  {
    char* ptr = (char*) VirtualAlloc(NULL,bytes,MEM_RESERVE,PAGE_READWRITE);
    if (ptr == NULL) throw std::runtime_error("memory reserve failed");
    return ptr;
  }

  void os_commit (void* ptr, size_t bytes) {
    VirtualAlloc(ptr,bytes,MEM_COMMIT,PAGE_READWRITE);
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

#if defined(__MIC__)

static double getFrequencyInMHz()
{
  struct timeval tvstart, tvstop;
  unsigned long long int cycles[2];
    
  gettimeofday(&tvstart, NULL);
  cycles[0] = rdtsc();
  gettimeofday(&tvstart, NULL);
  usleep(250000);
  gettimeofday(&tvstop, NULL);
  cycles[1] = rdtsc();
  gettimeofday(&tvstop, NULL);
  
  const unsigned long microseconds = ((tvstop.tv_sec-tvstart.tv_sec)*1000000) + (tvstop.tv_usec-tvstart.tv_usec);
  unsigned long mhz = (unsigned long) (cycles[1]-cycles[0]) / microseconds;

  //PRINT(mhz);
  return (double)mhz;
}

  static double micFrequency = getFrequencyInMHz();

#endif

namespace embree
{
void* os_malloc(size_t bytes)
  {
    int flags = MAP_PRIVATE | MAP_ANON;
#if defined(__MIC__)
    if (bytes > 16*4096) {
      flags |= MAP_HUGETLB | MAP_POPULATE;
      bytes = (bytes+2*1024*1024-1)&(-2*1024*1024);
    } else {
      bytes = (bytes+4095)&(-4096);
    }
#endif
    char* ptr = (char*) mmap(0, bytes, PROT_READ | PROT_WRITE, flags, -1, 0);
    if (ptr == NULL || ptr == MAP_FAILED) throw std::runtime_error("memory allocation failed");
    return ptr;
  }

  void* os_reserve(size_t bytes)
  {
    int flags = MAP_PRIVATE | MAP_ANON | MAP_NORESERVE;
#if defined(__MIC__)
    if (bytes > 16*4096) {
      flags |= MAP_HUGETLB;
      bytes = (bytes+2*1024*1024-1)&(-2*1024*1024);
    } else {
      bytes = (bytes+4095)&(-4096);
    }
#endif
    char* ptr = (char*) mmap(0, bytes, PROT_READ | PROT_WRITE, flags, -1, 0);
    if (ptr == NULL || ptr == MAP_FAILED) throw std::runtime_error("memory reserve failed");
    return ptr;
  }

  void os_commit (void* ptr, size_t bytes) {
  }

  void os_free(void* ptr, size_t bytes) 
  {
#if defined(__MIC__)
    if (bytes > 16*4096) {
      bytes = (bytes+2*1024*1024-1)&(-2*1024*1024);
    } else {
      bytes = (bytes+4095)&(-4096);
    }
#endif
    if (munmap(ptr,bytes) == -1) {
      throw std::runtime_error("memory free failed");
    }
  }


  double getSeconds() {
#if !defined(__MIC__)
    struct timeval tp; gettimeofday(&tp,NULL);
    return double(tp.tv_sec) + double(tp.tv_usec)/1E6;
#else
    return double(rdtsc()) / double(micFrequency*1E6);
#endif
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
