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

#ifndef __EMBREE_BVHI_BUILDER_FAST_H__
#define __EMBREE_BVHI_BUILDER_FAST_H__

#include "bvh4i.h"
#include "bvh4i_builder_util.h"
#include "bvh4i_builder_binner.h"

namespace embree
{
  namespace isa
  {
    class BVH4iBuilderFast : public Builder
    {
      ALIGNED_CLASS;
      
    public:
      
      /*! Constructor. */
      BVH4iBuilderFast (BVH4i* bvh, BuildSource* source, void* geometry, const size_t minLeafSize = 1, const size_t maxLeafSize = inf);
      
      /* build function */
      void build(size_t threadIndex, size_t threadCount);
      
      void allocateData();
      void freeData();
      
    public:
      TASK_FUNCTION(BVH4iBuilderFast,computePrimRefs);
      TASK_FUNCTION(BVH4iBuilderFast,buildSubTrees);
      TASK_FUNCTION(BVH4iBuilderFast,createTriangle1);
      TASK_FUNCTION(BVH4iBuilderFast,convertToSOALayout);
      TASK_RUN_FUNCTION(BVH4iBuilderFast,build_parallel);
      
    public:
      
      /*! build mode */
      enum { RECURSE = 1, BUILD_TOP_LEVEL = 3 };
      
      /*! splitting function that selects between sequential and parallel mode */
      bool split(BuildRecord& current, BuildRecord& left, BuildRecord& right, const size_t mode, const size_t threadID, const size_t numThreads);
      
      /*! perform sequential binning and splitting */
      bool splitSequential(BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild);
      
      /*! perform parallel binning and splitting */
      bool splitParallel(BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild, const size_t threadID, const size_t threads);
      
      /*! creates a leaf node */
      void createLeaf(BuildRecord& current, size_t threadIndex, size_t threadCount);
      
      /*! select between recursion and stack operations */
      void recurse(BuildRecord& current, const size_t mode, const size_t threadID, const size_t numThreads);
      
      /*! recursive build function */
      void recurseSAH(BuildRecord& current, const size_t mode, const size_t threadID, const size_t numThreads);
      
    protected:
      BuildSource* source;          //!< input geometry
      void* geometry;               //!< input geometry
      BVH4i* bvh;                   //!< Output BVH
      const PrimitiveType& primTy;  //!< triangle type stored in BVH
      
      /* work record handling */
    protected:
      static const size_t SIZE_WORK_STACK = 64;
      WorkStack<BuildRecord,SIZE_WORK_STACK> global_workStack;
      WorkStack<BuildRecord,SIZE_WORK_STACK> thread_workStack[MAX_MIC_THREADS];
      TaskScheduler::Task task;
      LinearBarrierActive global_barrier;
      ParallelBinner<16> parallelBinner;  
      
    protected:
      PrimRef* prims;
      BVHNode* node;
      Triangle1* accel;
      
    protected:
      size_t numPrimitives;
      size_t numNodes;
      size_t numAllocatedNodes;
      Centroid_Scene_AABB global_bounds;
      
      /*! node allocator */
    protected:
      AlignedAtomicCounter32  atomicID;
      
      __forceinline unsigned int allocNode(int size)
      {
        const unsigned int currentIndex = atomicID.add(size);
        if (unlikely(currentIndex >= numAllocatedNodes)) {
          FATAL("not enough nodes allocated");
        }
        return currentIndex;
      }
    };
  }
}

#endif
  
