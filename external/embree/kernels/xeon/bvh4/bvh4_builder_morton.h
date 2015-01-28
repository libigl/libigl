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

#ifndef __EMBREE_BVH4_BUILDER_MORTON_H__
#define __EMBREE_BVH4_BUILDER_MORTON_H__

#include "bvh4.h"
#include "../bvh4i/bvh4i_builder_util.h"

namespace embree
{
  namespace isa
  {
    class BVH4BuilderMorton : public Builder
    {
      ALIGNED_CLASS;
      
      /*! Type shortcuts */
      typedef BVH4::Node    Node;
      typedef BVH4::NodeRef NodeRef;
      typedef GlobalAllocator::ThreadAllocator Allocator;
      
      enum { RECURSE = 1, CREATE_TOP_LEVEL = 2 };
      
      static const size_t MAX_TOP_LEVEL_BINS = 1024;
      static const size_t NUM_TOP_LEVEL_BINS = 1024 + 4*BVH4::maxBuildDepth;

      static const size_t MORTON_LEAF_THRESHOLD = 4;
      static const size_t LATTICE_BITS_PER_DIM = 10;
      static const size_t LATTICE_SIZE_PER_DIM = size_t(1) << LATTICE_BITS_PER_DIM;
      
      static const size_t RADIX_BITS = 11;
      static const size_t RADIX_BUCKETS = (1 << RADIX_BITS);
      static const size_t RADIX_BUCKETS_MASK = (RADIX_BUCKETS-1);

    public:
  
      class __aligned(16) SmallBuildRecord 
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
        
        __forceinline bool operator<(const SmallBuildRecord& br) const { return size() < br.size(); } 
        __forceinline bool operator>(const SmallBuildRecord& br) const { return size() > br.size(); } 
      };

      struct __aligned(8) MortonID32Bit
      {
        union {
          struct {
            unsigned int code;
            unsigned int index;
          };
          int64 all;
        };
        
        __forceinline unsigned int get(const unsigned int shift, const unsigned and_mask) const {
          return (code >> shift) & and_mask;
        }
        
        __forceinline void operator=(const MortonID32Bit& v) {   
          all = v.all; 
        };  
        
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
          numBuildRecords = 0;
          numThreads = getNumberOfLogicalThreads();
          startGroup = new unsigned int[numThreads];
          startGroupOffset = new unsigned int[numThreads];
          radixCount = (ThreadRadixCountTy*) alignedMalloc(numThreads*sizeof(ThreadRadixCountTy));
        }

        ~MortonBuilderState () 
        {
          alignedFree(radixCount);
          delete[] startGroupOffset;
          delete[] startGroup;
        }

        size_t numThreads;
        unsigned int* startGroup;
        unsigned int* startGroupOffset;
        ThreadRadixCountTy* radixCount;
        
        size_t numBuildRecords;
        __aligned(64) SmallBuildRecord buildRecords[NUM_TOP_LEVEL_BINS];
        __aligned(64) WorkStack<SmallBuildRecord,NUM_TOP_LEVEL_BINS> workStack;
        LinearBarrierActive barrier;
      };
      
      /*! Constructor. */
      BVH4BuilderMorton (BVH4* bvh, BuildSource* source, Scene* scene, TriangleMeshScene::TriangleMesh* mesh, const size_t minLeafSize = 1, const size_t maxLeafSize = inf);
      
      /*! Destruction */
      ~BVH4BuilderMorton ();
      
      /* build function */
      void build(size_t threadIndex, size_t threadCount);
      
      /*! initialized the builder */
      void init(size_t threadIndex, size_t threadCount);
      
      /*! precalculate some per thread data */
      void initThreadState(const size_t threadID, const size_t numThreads);
      
      /*! single threaded build */
      void build_sequential_morton(size_t threadIndex, size_t threadCount);

      Centroid_Scene_AABB computeBounds();

      void computeMortonCodes(const size_t startID, const size_t endID, 
                              const size_t startGroup, const size_t startOffset, 
                              MortonID32Bit* __restrict__ const dest);

      /*! main build task */
      TASK_RUN_FUNCTION(BVH4BuilderMorton,build_parallel_morton);
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
      
      typedef void (*createLeafFunction)(const BVH4BuilderMorton* This, SmallBuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3f& box_o);
      typedef BBox3f (*leafBoundsFunction)(NodeRef& ref);
      
      /*! creates a leaf node */
      static void createTriangle1Leaf(const BVH4BuilderMorton* This, SmallBuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3f& box_o);
      static void createTriangle4Leaf(const BVH4BuilderMorton* This, SmallBuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3f& box_o);
      static void createTriangle1vLeaf(const BVH4BuilderMorton* This, SmallBuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3f& box_o);
      static void createTriangle4vLeaf(const BVH4BuilderMorton* This, SmallBuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3f& box_o);
      
      BBox3f createLeaf(SmallBuildRecord& current, Allocator& nodeAlloc, Allocator& leafAlloc, size_t threadID);
      
      /*! fallback split mode */
      void split_fallback(SmallBuildRecord& current, SmallBuildRecord& leftChild, SmallBuildRecord& rightChild) const;
      
      /*! split a build record into two */
      void split(SmallBuildRecord& current, SmallBuildRecord& left, SmallBuildRecord& right) const;
      
      /*! main recursive build function */
      BBox3f recurse(SmallBuildRecord& current, 
                     Allocator& nodeAlloc, Allocator& leafAlloc,
                     const size_t mode, 
                     const size_t threadID);
      
      static BBox3f leafBoundsTriangle1(NodeRef& ref);
      static BBox3f leafBoundsTriangle4(NodeRef& ref);
      static BBox3f leafBoundsTriangle1v(NodeRef& ref);
      static BBox3f leafBoundsTriangle4v(NodeRef& ref);
      BBox3f node_bounds(NodeRef& ref) const;
      
      /*! refit the toplevel part of the BVH */
      BBox3f refit_toplevel(NodeRef& index) const;
      
      /*! refit the sub-BVHs */
      BBox3f refit(NodeRef& index) const;
      
      /*! recreates morton codes when reaching a region where all codes are identical */
      void recreateMortonCodes(SmallBuildRecord& current) const;
      
    public:
      BVH4* bvh;               //!< Output BVH
      BuildSource* source;      //!< input geometry
      Scene* scene;
      TriangleMeshScene::TriangleMesh* mesh;
      
      size_t topLevelItemThreshold;
      size_t encodeShift;
      size_t encodeMask;

      static std::auto_ptr<MortonBuilderState> g_state;
            
    protected:
      MortonID32Bit* __restrict__ morton;
      size_t bytesMorton;
      
    public:
      size_t numGroups;
      size_t numPrimitives;
      size_t numAllocatedPrimitives;
      size_t numAllocatedNodes;
      Centroid_Scene_AABB global_bounds;
      createLeafFunction createSmallLeaf;
      leafBoundsFunction leafBounds;
      LockStepTaskScheduler scheduler;
      
    protected:
      __aligned(64) GlobalAllocator nodeAllocator;
      __aligned(64) GlobalAllocator primAllocator;
    };
  }
}
#endif
  
