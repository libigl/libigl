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

#include "sys/sysinfo.h"
#include "sys/sync/mutex.h"
#include "sys/taskscheduler.h"
#include "math/math.h"

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
    //enum { blockSize = 512*4096 };
    enum { blockSize = 16*4096 };
    //enum { blockSize = 4*4096 };
    
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
      clear();
    }

    /*! clears the allocator */
    void clear () 
    {
      for (size_t i=0; i<blocks.size(); i++) {
        Alloc::global.free(blocks[i]); 
      }
      ptr = NULL;
      cur = end = 0;
      blocks.resize(0);
    }

    /*! returns number of bytes allocated */
    size_t bytes () {
      return blocks.size() * Alloc::blockSize;
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
  class LinearAllocatorPerThread : public RefCount
  {
    ALIGNED_CLASS;

  public:

     /*! each thread handles block of that many bytes locally */
    enum { allocBlockSize = 4096 };

    /*! Per thread structure holding the current memory block. */
    struct __aligned(64) ThreadAllocator 
    {
      ALIGNED_CLASS_(64);
    public:

       /*! each thread handles block of that many bytes locally */
      enum { blockSize = allocBlockSize };

      /*! Default constructor. */
      __forceinline ThreadAllocator (LinearAllocatorPerThread* alloc = NULL) 
	: alloc(alloc), ptr(NULL), cur(0), end(0) {}

      /* Allocate aligned memory from the threads memory block. */
      __forceinline void* malloc(size_t bytes, size_t align = 16) 
      {
        cur += bytes + ((align - cur) & (align-1));
        if (likely(cur <= end)) return &ptr[cur - bytes];
        ptr = (char*) alloc->block.malloc(allocBlockSize);
        cur = 0;
        end = allocBlockSize;
        if (bytes > allocBlockSize) 
          THROW_RUNTIME_ERROR("allocated block is too large");
        cur += bytes;
        return &ptr[cur - bytes];
      }

      /*! clears the allocator */
      void clear () {
        ptr = NULL;
        cur = end = 0;
      }

    public:
      LinearAllocatorPerThread* alloc;
      char*  ptr;      //!< pointer to memory block
      size_t cur;      //!< Current location of the allocator.
      size_t end;      //!< End of the memory block.
    };

    /*! Allocator default construction. */
    LinearAllocatorPerThread () {}

    /*! Return pointer to start of memory region */
    __forceinline       void* base()       { return block.ptr; }
    __forceinline const void* base() const { return block.ptr; }
    __forceinline       void* curPtr()     { return block.ptr+block.cur; }

    /*! clears the allocator */
    void clear () {
      block.clear();
    }

    /*! initializes the allocator */
    void init (size_t bytesAllocate, size_t bytesReserve) 
    {
      clear();
      const size_t numThreads = getNumberOfLogicalThreads(); // FIXME: should get passed from outside
      bytesReserve = max(bytesAllocate,bytesReserve);
      size_t bytesReserved = max(bytesReserve,size_t(allocBlockSize*numThreads));
      block.init(bytesAllocate,bytesReserved);
    }

    /*! returns number of committed bytes */
    size_t bytes () const {
      return block.cur;
    }

    void shrink () {
      block.shrink();
    }

    void print_statistics()
    {
      size_t bytesAllocated = block.getAllocatedBytes();
      size_t bytesReserved = block.getReservedBytes();
      size_t bytesUsed = block.getUsedBytes();
      size_t bytesFree = block.getFreeBytes();
      
      printf("allocated = %3.2fMB, reserved = %3.2fMB, used = %3.2fMB (%3.2f%%), free = %3.2fMB (%3.2f%%)\n",
	     1E-6f*bytesAllocated, 1E-6f*bytesReserved,
	     1E-6f*bytesUsed, 100.0f*bytesUsed/bytesAllocated,
	     1E-6f*bytesFree, 100.0f*bytesFree/bytesAllocated);
    }

  private:

    struct Block 
    {
      Block () 
      : ptr(NULL), cur(0), reserveEnd(0), allocEnd(0), next(NULL) {}
      
      Block (size_t bytes, Block* next = NULL) 
      : ptr(NULL), cur(0), reserveEnd(bytes), allocEnd(0), next(next) {}

      ~Block () {
	if (ptr) os_free(ptr,reserveEnd); ptr = NULL;
	cur = reserveEnd = 0;
	if (next) delete next; next = NULL;
      }

      __forceinline void init (size_t bytesAllocate, size_t bytesReserved)
      {
	if (bytesReserved != size_t(reserveEnd) || bytesAllocate != allocEnd) 
	{
	  allocEnd = bytesAllocate;
	  if (ptr) os_free(ptr,reserveEnd);
	  ptr = (char*) os_reserve(bytesReserved);
	  os_commit(ptr,allocEnd);
	  reserveEnd = bytesReserved;
	}
      }

      __forceinline void clear() {
	cur = 0;
      }

      /*! Allocates some number of bytes. */
      void* malloc(size_t bytes) 
      {
	ssize_t i = atomic_add(&cur,bytes);
	if (unlikely(i+(ssize_t)bytes > reserveEnd)) THROW_RUNTIME_ERROR("build out of memory");
	void* p = &ptr[i];
	if (i+(ssize_t)bytes > allocEnd)
	  os_commit(p,bytes);
	return p;
      }

      void shrink () {
	if (ptr == NULL) return;
	os_shrink(ptr,cur,reserveEnd);
	reserveEnd = cur;
	allocEnd = cur;
      }

      size_t getAllocatedBytes() const {
	return allocEnd;
      }

      size_t getReservedBytes() const {
	return reserveEnd;
      }

      size_t getUsedBytes() const {
	return cur;
      }

      size_t getFreeBytes() const {
	return allocEnd-cur;
      }

    public:
      char*  ptr;                //!< pointer to memory
      atomic_t cur;              //!< Current location of the allocator.
      atomic_t allocEnd;
      atomic_t reserveEnd;              //!< End of the memory block.
      Block* next;
    };

  private:
    Block block;
  };



  class FastAllocator 
  {
    /*! maximal supported alignment */
    static const size_t maxAlignment = 64;

    /*! maximal allocation size */
    static const size_t maxAllocationSize = 2*1024*1024-maxAlignment;

  public:

    /*! Per thread structure holding the current memory block. */
    struct __aligned(64) Thread 
    {
      ALIGNED_CLASS_(64);
    public:

      /*! Constructor for usage with ThreadLocal */
      __forceinline Thread (void* alloc) 
	: alloc((FastAllocator*)alloc), ptr(NULL), cur(0), end(0), allocBlockSize(4096), bytesUsed(0), bytesWasted(0) {}

      /*! Default constructor. */
      __forceinline Thread (FastAllocator* alloc, const size_t allocBlockSize = 4096) 
	: alloc(alloc), ptr(NULL), cur(0), end(0), allocBlockSize(allocBlockSize), bytesUsed(0), bytesWasted(0)  {}

      /*! resets the allocator */
      __forceinline void reset() 
      {
	ptr = NULL;
	cur = end = 0;
	bytesWasted = bytesUsed = 0;
      }

      /* Allocate aligned memory from the threads memory block. */
      __forceinline void* malloc(size_t bytes, size_t align = 16) 
      {
        assert(align <= maxAlignment);
	bytesUsed += bytes;
	
        /* try to allocate in local block */
	size_t ofs = (align - cur) & (align-1); 
        cur += bytes + ofs;
        if (likely(cur <= end)) { bytesWasted += ofs; return &ptr[cur - bytes]; }
	cur -= bytes + ofs;

        /* if allocation is too large allocate with parent allocator */
        if (4*bytes > allocBlockSize) {
          return alloc->malloc(bytes,maxAlignment);
	}

#if 0 // FIXME: this optimization is broken

        /* get new partial block if allocation failed */
	if (alloc->usedBlocks) 
	{
	  size_t blockSize = allocBlockSize;
	  ptr = (char*) alloc->usedBlocks->malloc_some(blockSize,maxAlignment);
	  bytesWasted += end-cur;
	  cur = 0; end = blockSize;
	  
	  /* retry allocation */
	  size_t ofs = (align - cur) & (align-1); 
	  cur += bytes + ofs;
	  if (likely(cur <= end)) { bytesWasted += ofs; return &ptr[cur - bytes]; }
	  cur -= bytes + ofs;
	}
#endif

        /* get new full block if allocation failed */
        size_t blockSize = allocBlockSize;
	ptr = (char*) alloc->malloc(blockSize,maxAlignment);
	bytesWasted += end-cur;
	cur = 0; end = blockSize;
	
        /* retry allocation */
	ofs = (align - cur) & (align-1); 
        cur += bytes + ofs;
        if (likely(cur <= end)) { bytesWasted += ofs; return &ptr[cur - bytes]; }
	cur -= bytes + ofs;
	
        /* should never happen as large allocations get handled specially above */
        assert(false);
        return NULL;
      }

      /*! returns amount of used bytes */
      size_t getUsedBytes() const { return bytesUsed; }
      
      /*! returns amount of wasted bytes */
      size_t getWastedBytes() const { return bytesWasted + (end-cur); }

    public:
      FastAllocator* alloc;  //!< parent allocator
      char*  ptr;            //!< pointer to memory block
      size_t cur;            //!< current location of the allocator
      size_t end;            //!< end of the memory block
      size_t allocBlockSize; //!< block size for allocations
    private:
      size_t bytesWasted;    //!< number of bytes wasted
      size_t bytesUsed; //!< bumber of total bytes allocated
    };

    FastAllocator () 
      : growSize(4096), usedBlocks(NULL), freeBlocks(NULL), thread_local_allocators(this) {}

    ~FastAllocator () { 
      if (usedBlocks) usedBlocks->~Block(); usedBlocks = NULL;
      if (freeBlocks) freeBlocks->~Block(); freeBlocks = NULL;
    }

    /*! returns a fast thread local allocator */
    __forceinline Thread* instance() {
      return thread_local_allocators.get();
    }

    /*! initializes the allocator */
    void init(size_t bytesAllocate, size_t bytesReserve) {
      usedBlocks = Block::create(bytesAllocate,bytesReserve);
      growSize = bytesReserve;
    }

    /*! resets the allocator, memory blocks get reused */
    void reset () 
    {
      /* first reset all used blocks */
      if (usedBlocks) usedBlocks->reset();

      /* find end of free block list */
      Block* volatile& freeBlocksEnd = freeBlocks;
      while (freeBlocksEnd) freeBlocksEnd = freeBlocksEnd->next;

      /* add previously used blocks to end of free block list */
      freeBlocksEnd = usedBlocks;
      usedBlocks = NULL;

      /* reset all thread local allocators */
      thread_local_allocators.reset();
    }

    /*! shrinks all memory blocks to the actually used size */
    void shrink () {
      usedBlocks->shrink();
      if (freeBlocks) freeBlocks->~Block(); freeBlocks = NULL;
    }

    /*! thread safe allocation of memory */
    void* malloc(size_t bytes, size_t align) 
    {
      assert(align <= maxAlignment);

      while (true) 
      {
        /* allocate using current block */
	Block* myUsedBlocks = usedBlocks;
        if (myUsedBlocks) {
          void* ptr = usedBlocks->malloc(bytes,align);
          if (ptr) return ptr;
        }

        /* throw error if allocation is too large */
        if (bytes > maxAllocationSize)
          THROW_RUNTIME_ERROR("allocation is too large");

        /* if this fails allocate new block */
        {
          Lock<AtomicMutex> lock(mutex);
	  if (myUsedBlocks == usedBlocks)
	  {
	    if (freeBlocks) {
	      Block* nextFreeBlock = freeBlocks->next;
	      freeBlocks->next = usedBlocks;
	      __memory_barrier();
	      usedBlocks = freeBlocks;
	      freeBlocks = nextFreeBlock;
	    } else {
	      growSize = min(2*growSize,size_t(maxAllocationSize+maxAlignment));
	      usedBlocks = Block::create(growSize-maxAlignment, growSize-maxAlignment, usedBlocks);
	    }
	  }
        }
      }
    }

    void print_statistics()
    {
      size_t bytesFree = 0;
      size_t bytesAllocated = 0;
      size_t bytesReserved = 0;
      size_t bytesUsed = 0;
      size_t bytesWasted = 0;

      if (freeBlocks) {
	bytesFree += freeBlocks->getAllocatedBytes();
	bytesAllocated += freeBlocks->getAllocatedBytes();
	bytesReserved += freeBlocks->getReservedBytes();
      }
      if (usedBlocks) {
	bytesFree += usedBlocks->getFreeBytes();
	bytesAllocated += usedBlocks->getAllocatedBytes();
	bytesReserved += usedBlocks->getReservedBytes();
	
	Block* cur = usedBlocks;
	while ((cur = cur->next) != NULL)
	  bytesWasted += cur->getFreeBytes();
      }

      for (size_t t=0; t<thread_local_allocators.threads.size(); t++) {
	bytesUsed   += thread_local_allocators.threads[t]->getUsedBytes();
	bytesWasted += thread_local_allocators.threads[t]->getWastedBytes();
      }
      
      printf("allocated = %3.2fMB, reserved = %3.2fMB, used = %3.2fMB (%3.2f%%), wasted = %3.2fMB (%3.2f%%), free = %3.2fMB (%3.2f%%)\n",
	     1E-6f*bytesAllocated, 1E-6f*bytesReserved,
	     1E-6f*bytesUsed, 100.0f*bytesUsed/bytesAllocated,
	     1E-6f*bytesWasted, 100.0f*bytesWasted/bytesAllocated,
	     1E-6f*bytesFree, 100.0f*bytesFree/bytesAllocated);
    }

  private:

    struct Block 
    {
      static Block* create(size_t bytesAllocate, size_t bytesReserve, Block* next = NULL)
      {
        void* ptr = os_reserve(sizeof(Block)+bytesReserve);
        os_commit(ptr,sizeof(Block)+bytesAllocate);
        bytesAllocate = ((sizeof(Block)+bytesAllocate+4095) & ~(4095)) - sizeof(Block); // always comsume full pages
        bytesReserve  = ((sizeof(Block)+bytesReserve +4095) & ~(4095)) - sizeof(Block); // always comsume full pages
        return new (ptr) Block(bytesAllocate,bytesReserve,next);
      }

      Block (size_t bytesAllocate, size_t bytesReserve, Block* next) 
      : cur(0), allocEnd(bytesAllocate), reserveEnd(bytesReserve), next(next) {}

      ~Block () {
	if (next) next->~Block(); next = NULL;
        os_free(this,sizeof(Block)+reserveEnd);
      }

      void* malloc(size_t bytes, size_t align = 16) 
      {
        assert(align <= maxAlignment);
        bytes = (bytes+(align-1)) & ~(align-1); // FIXME: works only if all alignments are equal
	if (unlikely(cur+bytes > reserveEnd)) return NULL;
	const size_t i = atomic_add(&cur,bytes);
	if (unlikely(i+bytes > reserveEnd)) return NULL;
	if (i+bytes > allocEnd) os_commit(&data[i],bytes); // FIXME: optimize, may get called frequently
	return &data[i];
      }
      
      void* malloc_some(size_t& bytes, size_t align = 16) 
      {
        assert(align <= maxAlignment);
        bytes = (bytes+(align-1)) & ~(align-1); // FIXME: works only if all alignments are equal
	const size_t i = atomic_add(&cur,bytes);
	if (unlikely(i+bytes > reserveEnd)) bytes = reserveEnd-i;
	if (i+bytes > allocEnd) os_commit(&data[i],bytes); // FIXME: optimize, may get called frequently
	return &data[i];
      }

      void reset () 
      {
        allocEnd = max(allocEnd,(size_t)cur);
        cur = 0;
        if (next) next->reset();
      }

      void shrink () 
      {
        os_shrink(&data[0],cur,reserveEnd);
        reserveEnd = allocEnd = cur;
        if (next) next->shrink();
      }

      size_t getAllocatedBytes() const {
	return allocEnd + (next ? next->getAllocatedBytes() : 0);
      }

      size_t getReservedBytes() const {
	return reserveEnd + (next ? next->getReservedBytes() : 0);
      }

      size_t getFreeBytes() const {
	return allocEnd-cur;
      }

    public:
      atomic_t cur;              //!< current location of the allocator
      size_t allocEnd;           //!< end of the allocated memory region
      size_t reserveEnd;         //!< end of the reserved memory region
      Block* next;               //!< pointer to next block in list
      char align[maxAlignment-4*sizeof(size_t)]; //!< align data to maxAlignment
      char data[];               //!< here starts memory to use for allocations
    };

  private:
    AtomicMutex mutex;
    Block* volatile usedBlocks;
    Block* volatile freeBlocks;
    size_t growSize;

    ThreadLocal<Thread> thread_local_allocators; //!< thread local allocators

  private:
    size_t bytesWasted;    //!< number of bytes wasted
    size_t bytesUsed; //!< bumber of total bytes allocated
  };
}
