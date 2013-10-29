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

#ifndef __EMBREE_PRIMREF_ALLOCATOR_H__
#define __EMBREE_PRIMREF_ALLOCATOR_H__

#include "../common/alloc.h"
#include "primrefblock.h"

namespace embree
{
  class PrimRefAlloc : public AllocatorBase
  {
    ALIGNED_CLASS;
  public:
   
    struct __align(4096) ThreadPrimBlockAllocator 
    {
      ALIGNED_CLASS_(4096);
    public:
      
      __forceinline atomic_set<PrimRefBlock>::item* malloc(size_t thread, AllocatorBase* alloc) 
      {
        /* try to take a block from local list */
        atomic_set<PrimRefBlock>::item* ptr = local_free_blocks.take_unsafe();
        if (ptr) return new (ptr) atomic_set<PrimRefBlock>::item();
        
        /* if this failed again we have to allocate more memory */
        ptr = (atomic_set<PrimRefBlock>::item*) alloc->malloc(sizeof(atomic_set<PrimRefBlock>::item));
        
        /* return first block */
        return new (ptr) atomic_set<PrimRefBlock>::item();
      }
      
      __forceinline void free(atomic_set<PrimRefBlock>::item* ptr) {
        local_free_blocks.insert_unsafe(ptr);
      }

    public:
      atomic_set<PrimRefBlock> local_free_blocks; //!< only accessed from one thread
    };

  public:

    /*! Allocator default construction. */
    PrimRefAlloc () {
      threadPrimBlockAllocator = new ThreadPrimBlockAllocator[getNumberOfLogicalThreads()];
    }

    /*! Allocator destructor. */
    virtual ~PrimRefAlloc() {
      delete[] threadPrimBlockAllocator; threadPrimBlockAllocator = NULL;
    }

    /*! Allocate a primitive block */
    __forceinline atomic_set<PrimRefBlock>::item* malloc(size_t thread) {
      return threadPrimBlockAllocator[thread].malloc(thread,this);
    }

    /*! Frees a primitive block */
    __forceinline void free(size_t thread, atomic_set<PrimRefBlock>::item* block) {
      return threadPrimBlockAllocator[thread].free(block);
    }

  private:
    ThreadPrimBlockAllocator* threadPrimBlockAllocator;  //!< Thread local allocator
  };
}

#endif
