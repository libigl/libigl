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

#include "geometry/primitive.h"
#include "geometry/bezier1v.h"
#include "builders/primrefalloc.h"
#include "builders/heuristic_split.h"

namespace embree
{
  namespace isa
  {
    class BVH4BuilderHair : public Builder
    {
      ALIGNED_CLASS;
    public:
      
      /*! builder entry point */
      void build(size_t threadIndex, size_t threadCount);
      
      /*! Constructor. */
      BVH4BuilderHair (BVH4* bvh, Scene* scene, size_t mode);
      
    protected:
      
      typedef BezierPrim PrimRef;
      typedef PrimRefBlockT<BezierPrim> PrimRefBlock;
      typedef atomic_set<PrimRefBlock> BezierRefList;
      typedef LinearAllocatorPerThread::ThreadAllocator Allocator;
      
      /*! stores all info to build a subtree */
      struct BuildTask
      {
	__forceinline BuildTask () {}
	
	__forceinline BuildTask (BVH4::NodeRef* dst, size_t depth, BezierRefList& prims, const PrimInfo& pinfo, const NAABBox3fa& bounds, const PrimInfo& sinfo, const Split& split)
	  : dst(dst), depth(depth), prims(prims), pinfo(pinfo), bounds(bounds), sinfo(sinfo), split(split) {}

      public:
	__forceinline friend bool operator< (const BuildTask& a, const BuildTask& b) {
	  //return halfArea(a.bounds.bounds) < halfArea(b.bounds.bounds);
	  return a.pinfo.size() < b.pinfo.size();
	}
	
      public:
	BVH4::NodeRef* dst;
	size_t depth;
	BezierRefList prims;
	PrimInfo pinfo;
	NAABBox3fa bounds;
	PrimInfo sinfo;
	Split split;
      };
      
    private:
      
      /*! creates a leaf node */
      virtual BVH4::NodeRef createLeaf(size_t threadIndex, Allocator& nodeAlloc, Allocator& leafAlloc, size_t depth, BezierRefList& prims, const PrimInfo& pinfo) = 0;
      
      /*! creates a large leaf that could be larger than supported by the BVH */
      BVH4::NodeRef createLargeLeaf(size_t threadIndex, Allocator& nodeAlloc, Allocator& leafAlloc, BezierRefList& prims, const PrimInfo& pinfo, size_t depth);
    
      template<bool Parallel>
	Split find_split(size_t threadIndex, size_t threadCount, BezierRefList& prims, const PrimInfo& pinfo, const NAABBox3fa& bounds, const PrimInfo& sinfo);
      
      /*! execute single task and create subtasks */
      template<bool Parallel>
	void processTask(size_t threadIndex, size_t threadCount, Allocator& nodeAlloc, Allocator& leafAlloc, BuildTask& task, BuildTask task_o[BVH4::N], size_t& N);
      
      /*! recursive build function for aligned and non-aligned bounds */
      void recurseTask(size_t threadIndex, size_t threadCount, Allocator& nodeAlloc, Allocator& leafAlloc, BuildTask& task);
      
      TASK_SET_FUNCTION(BVH4BuilderHair,task_build_parallel);
      
    public:
      Scene* scene;          //!< source
      size_t minLeafSize;    //!< minimal size of a leaf
      size_t maxLeafSize;    //!< maximal size of a leaf
      bool enableSpatialSplits; //!< turns on spatial splits
      size_t listMode;
      
      BVH4* bvh;         //!< output
      LockStepTaskScheduler* scheduler;
      PrimRefBlockAlloc<PrimRef> alloc;                 //!< Allocator for primitive blocks
      
      MutexSys taskMutex;
      volatile atomic_t numActiveTasks;
      volatile atomic_t numGeneratedPrims;
      volatile atomic_t remainingReplications;
      vector_t<BuildTask> tasks;
    };

    /*! specializes the builder for different leaf types */
    template<typename Triangle>
    class BVH4BuilderHairT : public BVH4BuilderHair
    {
    public:
      BVH4BuilderHairT (BVH4* bvh, Scene* scene, size_t mode);
      BVH4::NodeRef createLeaf(size_t threadIndex, Allocator& nodeAlloc, Allocator& leafAlloc, size_t depth, BezierRefList& prims, const PrimInfo& pinfo);
    };
  }
}
