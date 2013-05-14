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

#ifndef __EMBREE_ALLOC_H__
#define __EMBREE_ALLOC_H__

#include "sys/sysinfo.h"
#include "sys/sync/mutex.h"
#include "sys/taskscheduler.h"

#include <vector>

namespace embree
{
  /*! Global memory pool. Node, triangle, and intermediary build data
      is allocated from this memory pool and returned to it. The pool
      does not return memory to the operating system unless the clear function
      is called. */
  class Alloc
  {
  public:

    /*! Allocation block size. */
    enum { blockSize = 512*4096 };
    
    /*! single allocator object */
    static Alloc global;

    /*! Allocator default construction. */
    Alloc ();

    /*! Allocator destructor. */
    ~Alloc ();

    /*! returns size of memory pool */
    size_t size() const;
    
    /*! frees all available memory */
    void clear();
    
    /*! allocates a memory block */
    void* malloc();
    
    /*! frees a memory block */
    void free(void* ptr);
    
  private:
    MutexSys mutex;                 //<! Mutex to protect access to blocks vector
    std::vector<void*> blocks;      //<! list of available memory blocks
  };

  /*! Base class for a each memory allocator. Allocates from blocks of the 
    Alloc class and returns these blocks on destruction. */
  class AllocatorBase 
  {
  public:

    /*! Default constructor. */
    AllocatorBase () : ptr(NULL), cur(0), end(0) {
    }
    
    /*! Returns all allocated blocks to Alloc class. */
    ~AllocatorBase () {
      for (size_t i=0; i<blocks.size(); i++) {
        Alloc::global.free(blocks[i]); 
      }
    }

    /*! Allocates some number of bytes. */
    void* malloc(size_t bytes) 
    {
      Lock<MutexSys> lock(mutex);
      cur += bytes;
      if (cur <= end) return &ptr[cur - bytes];
      ptr = (char*) Alloc::global.malloc();
      blocks.push_back(ptr);
      cur = 0;
      end = Alloc::blockSize;
      assert(bytes<=Alloc::blockSize);
      cur += bytes;
      return &ptr[cur - bytes];
    }
    
  private:
    MutexSys mutex;                  //!< mutex to protect access to this class
    char*  ptr;                      //!< pointer to memory block
    size_t cur;                      //!< Current location of the allocator.
    size_t end;                      //!< End of the memory block.
    std::vector<void*> blocks;       //!< available memory blocks
  };

  /*! This class implements an efficient multi-threaded memory
   *  allocation scheme. The per thread allocator allocates from its
   *  current memory block or requests a new block from the slower
   *  global allocator when its block is full. */
  class AllocatorPerThread : public AllocatorBase
  {
    ALIGNED_CLASS;

     /*! Allocation block size. Number of bytes to request from the
      *  base allocator when the memory block of a thread is empty. */
    enum { allocBlockSize = 4096*16 };

  public:

    /*! Allocator default construction. */
    AllocatorPerThread () {
      thread = new ThreadAllocator[getNumberOfLogicalThreads()];
    }

    /*! Allocator destructor. */
    ~AllocatorPerThread() {
      delete[] thread; thread = NULL;
    }

    /*! Aligned memory allocation */
    __forceinline void* malloc(const TaskScheduler::ThreadInfo& tinfo, size_t bytes, size_t align = 16) {
      return thread[tinfo.id].malloc(bytes,align,this);
    }

  private:

     /*! Per thread structure holding the current memory block. */
    struct __align(4096) ThreadAllocator 
    {
      ALIGNED_CLASS_(4096);
    public:

      /*! Default constructor. */
      __forceinline ThreadAllocator () : ptr(NULL), cur(0), end(0) {}

      /* Allocate aligned memory from the threads memory block. */
      __forceinline void* malloc(size_t bytes, size_t align, AllocatorBase* alloc) 
      {
        cur += (align - cur) & (align-1);
        cur += bytes;
        if (cur <= end) return &ptr[cur - bytes];
        ptr = (char*) alloc->malloc(allocBlockSize);
        cur = 0;
        end = allocBlockSize;
        assert(bytes<=allocBlockSize);
        cur += bytes;
        return &ptr[cur - bytes];
      }

    public:
      char*  ptr;      //!< pointer to memory block
      size_t cur;      //!< Current location of the allocator.
      size_t end;      //!< End of the memory block.
    };

  private:
    ThreadAllocator* thread;   //!< one allocator for each thread
  };
}

#endif
