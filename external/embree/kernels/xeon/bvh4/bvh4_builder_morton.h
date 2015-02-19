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

#include "bvh4.h"
#include "builders/heuristic_fallback.h"
#include "builders/workstack.h"

namespace embree
{
  namespace isa
  {
    class BVH4BuilderMorton : public Builder
    {
      ALIGNED_CLASS;
      
    protected:
      /*! Type shortcuts */
      typedef BVH4::Node    Node;
      typedef BVH4::NodeRef NodeRef;
      //typedef GlobalAllocator::ThreadAllocator Allocator;
      typedef LinearAllocatorPerThread::ThreadAllocator Allocator;
      
      enum { RECURSE = 1, CREATE_TOP_LEVEL = 2 };
      
      static const size_t MAX_TOP_LEVEL_BINS = 1024;
      static const size_t NUM_TOP_LEVEL_BINS = 1024 + 4*BVH4::maxBuildDepth;

      static const size_t LATTICE_BITS_PER_DIM = 10;
      static const size_t LATTICE_SIZE_PER_DIM = size_t(1) << LATTICE_BITS_PER_DIM;
      
      static const size_t RADIX_BITS = 11;
      static const size_t RADIX_BUCKETS = (1 << RADIX_BITS);
      static const size_t RADIX_BUCKETS_MASK = (RADIX_BUCKETS-1);

    public:
  
      class BuildRecord 
      {
      public:
        unsigned int begin;
        unsigned int end;
        unsigned int depth;
        BVH4::NodeRef* parent;
        
        __forceinline unsigned int size() const {
          return end - begin;
        }
        
        __forceinline void init(const unsigned int _begin, const unsigned int _end)			 
        {
          begin = _begin;
          end = _end;
          depth = 1;
          parent = NULL;
	}
        
	struct Greater {
	  __forceinline bool operator()(const BuildRecord& a, const BuildRecord& b) {
	    return a.size() > b.size();
	  }
	};
      };

      struct __aligned(8) MortonID32Bit
      {
        union {
          struct {
	    unsigned int code;
	    unsigned int index;
	    //uint64 index;
          };
          //int64 all;
        };

        __forceinline operator unsigned() const { return code; }
        
        __forceinline unsigned int get(const unsigned int shift, const unsigned and_mask) const {
          return (code >> shift) & and_mask;
        }
        
        /*__forceinline void operator=(const MortonID32Bit& v) {   
          all = v.all; 
	  };*/  
        
        __forceinline friend std::ostream &operator<<(std::ostream &o, const MortonID32Bit& mc) {
          o << "index " << mc.index << " code = " << mc.code;
          return o;
        }
        
        __forceinline bool operator<(const MortonID32Bit &m) const { return code < m.code; } 
        __forceinline bool operator>(const MortonID32Bit &m) const { return code > m.code; } 
      };

      struct MortonBuilderState
      {
        ALIGNED_CLASS;

      public:

        typedef unsigned int ThreadRadixCountTy[RADIX_BUCKETS];

        MortonBuilderState () 
        {
	  taskCounter = 0;
          numThreads = getNumberOfLogicalThreads();
          startGroup = new unsigned int[numThreads];
          startGroupOffset = new unsigned int[numThreads];
	  dest = new size_t[numThreads];
          radixCount = (ThreadRadixCountTy*) alignedMalloc(numThreads*sizeof(ThreadRadixCountTy));
        }

        ~MortonBuilderState () 
        {
          alignedFree(radixCount);
          delete[] startGroupOffset;
          delete[] startGroup;
	  delete[] dest;
        }

        size_t numThreads;
        unsigned int* startGroup;
        unsigned int* startGroupOffset;
	size_t* dest;
        ThreadRadixCountTy* radixCount;
        
	atomic_t taskCounter;
	std::vector<BuildRecord> buildRecords;
        __aligned(64) WorkStack<BuildRecord,NUM_TOP_LEVEL_BINS> workStack;
        LinearBarrierActive barrier;
      };
      
      /*! Constructor. */
      BVH4BuilderMorton (BVH4* bvh, Scene* scene, TriangleMesh* mesh, size_t listMode, size_t logBlockSize, bool needVertices, size_t primBytes, const size_t minLeafSize, const size_t maxLeafSize);
      
      /*! Destruction */
      ~BVH4BuilderMorton ();
      
      /* build function */
      void build(size_t threadIndex, size_t threadCount);
      
      /*! precalculate some per thread data */
      void initThreadState(const size_t threadID, const size_t numThreads);
      
      /*! single threaded build */
      void build_sequential_morton(size_t threadIndex, size_t threadCount);

      CentGeomBBox3fa computeBounds();

      void computeMortonCodes(const size_t startID, const size_t endID, size_t& destID,
                              const size_t startGroup, const size_t startOffset, 
                              MortonID32Bit* __restrict__ const dest);

      /*! main build task */
      TASK_SET_FUNCTION(BVH4BuilderMorton,build_parallel_morton);
      TaskScheduler::Task task;
      
      /*! task that calculates the bounding box of the scene */
      TASK_FUNCTION(BVH4BuilderMorton,computeBounds);
      
      /*! task that calculates the morton codes for each primitive in the scene */
      TASK_FUNCTION(BVH4BuilderMorton,computeMortonCodes);
      
      /*! parallel sort of the morton codes */
      TASK_FUNCTION(BVH4BuilderMorton,radixsort);
      
      /*! task that builds a list of sub-trees */
      TASK_FUNCTION(BVH4BuilderMorton,recurseSubMortonTrees);
      
    public:
      
      /*! creates leaf node */
      virtual void createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3fa& box_o) = 0;
      
      /*! creates leaf node that is larger than supported by BVH */
      BBox3fa createLeaf(BuildRecord& current, Allocator& nodeAlloc, Allocator& leafAlloc, size_t threadID);
      
      /*! fallback split mode */
      void splitFallback(BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild) const;
      
      /*! split a build record into two */
      void split(BuildRecord& current, BuildRecord& left, BuildRecord& right) const;
      
      /*! main recursive build function */
      BBox3fa recurse(BuildRecord& current, 
                     Allocator& nodeAlloc, Allocator& leafAlloc,
                     const size_t mode, 
                     const size_t threadID);

      /*! calculates bounding box of leaf node */
      virtual BBox3fa leafBounds(NodeRef& ref) const = 0;

      /*! calculates bounding box of node */
      BBox3fa nodeBounds(NodeRef& ref) const;
      
      /*! refit the toplevel part of the BVH */
      BBox3fa refitTopLevel(NodeRef& index) const;
      
      /*! refit the sub-BVHs */
      BBox3fa refit(NodeRef& index) const;
      
      /*! recreates morton codes when reaching a region where all codes are identical */
      void recreateMortonCodes(BuildRecord& current) const;
      
    public:
      BVH4* bvh;               //!< Output BVH
      LockStepTaskScheduler* scheduler;
      std::unique_ptr<MortonBuilderState> state;

      Scene* scene;
      TriangleMesh* mesh;
      size_t logBlockSize;
      size_t blocks(size_t N) { return (N+((1<<logBlockSize)-1)) >> logBlockSize; }
      bool needVertices;
      size_t primBytes; 
      size_t minLeafSize;
      size_t maxLeafSize;
      size_t listMode;

      size_t topLevelItemThreshold;
      size_t encodeShift;
      size_t encodeMask;
            
    public:
      MortonID32Bit* __restrict__ morton;
      size_t bytesMorton;
      
    public:
      size_t numGroups;
      size_t numPrimitives;
      size_t numAllocatedPrimitives;
      size_t numAllocatedNodes;
      CentGeomBBox3fa global_bounds;
      Barrier barrier;
      //createSmallLeaf createSmallLeaf;
      //leafBounds leafBounds;
      //LockStepTaskScheduler scheduler;
    };

    class BVH4Triangle1BuilderMorton : public BVH4BuilderMorton
    {
    public:
      BVH4Triangle1BuilderMorton (BVH4* bvh, Scene* scene, size_t listMode);
      BVH4Triangle1BuilderMorton (BVH4* bvh, TriangleMesh* mesh, size_t listMode);
      BBox3fa leafBounds(NodeRef& ref) const;
      void createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3fa& box_o);
    };

    class BVH4Triangle4BuilderMorton : public BVH4BuilderMorton
    {
    public:
      BVH4Triangle4BuilderMorton (BVH4* bvh, Scene* scene, size_t listMode);
      BVH4Triangle4BuilderMorton (BVH4* bvh, TriangleMesh* mesh, size_t listMode);
      BBox3fa leafBounds(NodeRef& ref) const;
      void createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3fa& box_o);
    };

    class BVH4Triangle8BuilderMorton : public BVH4BuilderMorton
    {
    public:
      BVH4Triangle8BuilderMorton (BVH4* bvh, Scene* scene, size_t listMode);
      BVH4Triangle8BuilderMorton (BVH4* bvh, TriangleMesh* mesh, size_t listMode);
      BBox3fa leafBounds(NodeRef& ref) const;
      void createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3fa& box_o);
    };

    class BVH4Triangle1vBuilderMorton : public BVH4BuilderMorton
    {
    public:
      BVH4Triangle1vBuilderMorton (BVH4* bvh, Scene* scene, size_t listMode);
      BVH4Triangle1vBuilderMorton (BVH4* bvh, TriangleMesh* mesh, size_t listMode);
      BBox3fa leafBounds(NodeRef& ref) const;
      void createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3fa& box_o);
    };

    class BVH4Triangle4vBuilderMorton : public BVH4BuilderMorton
    {
    public:
      BVH4Triangle4vBuilderMorton (BVH4* bvh, Scene* scene, size_t listMode);
      BVH4Triangle4vBuilderMorton (BVH4* bvh, TriangleMesh* mesh, size_t listMode);
      BBox3fa leafBounds(NodeRef& ref) const;
      void createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3fa& box_o);
    };

    class BVH4Triangle4iBuilderMorton : public BVH4BuilderMorton
    {
    public:
      BVH4Triangle4iBuilderMorton (BVH4* bvh, Scene* scene, size_t listMode);
      BVH4Triangle4iBuilderMorton (BVH4* bvh, TriangleMesh* mesh, size_t listMode);
      BBox3fa leafBounds(NodeRef& ref) const;
      void createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3fa& box_o);
    };
  }
}
