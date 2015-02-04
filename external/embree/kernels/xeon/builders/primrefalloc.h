// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
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

#pragma once

#include "../common/alloc.h"
#include "primrefblock.h"

namespace embree
{
  template<typename PrimRef>
    class PrimRefBlockAlloc : public AllocatorBase
  {
    ALIGNED_CLASS;
  public:
    
    struct __aligned(4096) ThreadPrimBlockAllocator 
    {
      ALIGNED_CLASS_(4096);
    public:
      
      __forceinline typename atomic_set<PrimRefBlockT<PrimRef> >::item* malloc(size_t thread, AllocatorBase* alloc) 
      {
	/* try to take a block from local list */
	typename atomic_set<PrimRefBlockT<PrimRef> >::item* ptr = local_free_blocks.take_unsafe();
	if (ptr) return new (ptr) typename atomic_set<PrimRefBlockT<PrimRef> >::item();
	
	/* if this failed again we have to allocate more memory */
	//ptr = (typename atomic_set<PrimRefBlockT<PrimRef> >::item*) alloc->malloc_global(sizeof(typename atomic_set<PrimRefBlockT<PrimRef> >::item));
        ptr = (typename atomic_set<PrimRefBlockT<PrimRef> >::item*) alloc->malloc(sizeof(typename atomic_set<PrimRefBlockT<PrimRef> >::item));
        
	/* return first block */
	return new (ptr) typename atomic_set<PrimRefBlockT<PrimRef> >::item();
      }
      
      __forceinline void free(typename atomic_set<PrimRefBlockT<PrimRef> >::item* ptr) {
	local_free_blocks.insert_unsafe(ptr);
      }
      
    public:
      atomic_set<PrimRefBlockT<PrimRef> > local_free_blocks; //!< only accessed from one thread
    };
    
  public:
    
    /*! Allocator default construction. */
    PrimRefBlockAlloc () /*: ptr(NULL), cur(0), end(0), bytesAllocated(0)*/ {
      threadPrimBlockAllocator = new ThreadPrimBlockAllocator[getNumberOfLogicalThreads()];
    }
    
    /*! Allocator destructor. */
    virtual ~PrimRefBlockAlloc() {
      delete[] threadPrimBlockAllocator; threadPrimBlockAllocator = NULL;
      /*if (ptr) os_free(ptr,end); ptr = NULL;
	cur = end = 0;*/
    }
    
    /*! Allocate a primitive block */
    __forceinline typename atomic_set<PrimRefBlockT<PrimRef> >::item* malloc(size_t thread) {
      return threadPrimBlockAllocator[thread].malloc(thread,this);
    }
    
    /*! Frees a primitive block */
    __forceinline void free(size_t thread, typename atomic_set<PrimRefBlockT<PrimRef> >::item* block) {
      return threadPrimBlockAllocator[thread].free(block);
    }

#if 0
    /*! initializes the allocator */
    void init (size_t numAllocate, size_t numReserve) 
    {
      size_t bytesAllocate = numAllocate*sizeof(PrimRef);
      size_t bytesReserve = numReserve*sizeof(PrimRef);
      const size_t numThreads = getNumberOfLogicalThreads();
      bytesReserve = max(bytesAllocate,bytesReserve);
      size_t bytesReserved = max(bytesReserve,size_t(sizeof(typename atomic_set<PrimRefBlockT<PrimRef> >::item)*numThreads));
      if (bytesReserved != size_t(end) || bytesAllocate != bytesAllocated) 
      {
        bytesAllocated = bytesAllocate;
        if (ptr) os_free(ptr,end);
        ptr = (char*) os_reserve(bytesReserved);
        os_commit(ptr,bytesAllocated);
        //memset(ptr,0,bytesAllocated);
        end = bytesReserved;
      }
    }

  private:

    /*! Allocates some number of bytes. */
    void* malloc_global(size_t bytes) 
    {
      ssize_t i = atomic_add(&cur,bytes);
      if (unlikely(i > end)) THROW_RUNTIME_ERROR("build out of memory");
      void* p = &ptr[i];
      if (i+(ssize_t)bytes > bytesAllocated)
        os_commit(p,bytes);
      return p;
    }
#endif
    
  private:
    ThreadPrimBlockAllocator* threadPrimBlockAllocator;  //!< Thread local allocator

#if 0
    char*  ptr;                //!< pointer to memory
    atomic_t cur;              //!< Current location of the allocator.
    atomic_t end;              //!< End of the memory block.
    atomic_t bytesAllocated;
#endif
  };
  
  typedef PrimRefBlockAlloc<PrimRef> PrimRefAlloc;
}
