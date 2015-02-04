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
#include "builders/heuristic_object_partition.h"
#include "builders/workstack.h"
#include "common/scene_user_geometry.h"
#include "common/scene_subdiv_mesh.h"
#include "geometry/grid.h"

#include "algorithms/parallel_for_for.h"
#include "algorithms/parallel_for_for_prefix_sum.h"


namespace embree
{
  namespace isa
  {
    class BVH4BuilderFast : public Builder
    {
      ALIGNED_CLASS;

    protected:
      typedef BVH4::Node Node;
      typedef BVH4::NodeRef NodeRef;
      typedef LinearAllocatorPerThread::ThreadAllocator Allocator;
      static const size_t SIZE_WORK_STACK = 64;

    public:

      class __aligned(64) BuildRecord : public PrimInfo
      {
      public:
	unsigned int depth;         //!< depth from the root of the tree
	float sArea;
	BVH4::NodeRef* parent; 
	
        BuildRecord() {}

#if defined(_MSC_VER)
        BuildRecord& operator=(const BuildRecord &arg) { 
          memcpy(this, &arg, sizeof(BuildRecord));    
          return *this;
        }
#endif
	__forceinline void init(unsigned int depth)
	{
	  this->depth = depth;
	  sArea = area(geomBounds);
	}

	__forceinline void init(const CentGeomBBox3fa& _bounds, const unsigned int _begin, const unsigned int _end)
	{
	  geomBounds = _bounds.geomBounds;
	  centBounds = _bounds.centBounds;
	  begin  = _begin;
	  end    = _end;
	  sArea = area(geomBounds);
	}
	
	__forceinline float sceneArea() {
	  return sArea;
	}
	
	__forceinline bool operator<(const BuildRecord &br) const { return size() < br.size(); } 
	__forceinline bool operator>(const BuildRecord &br) const { return size() > br.size(); } 
	
	struct Greater {
	  bool operator()(const BuildRecord& a, const BuildRecord& b) {
	    return a > b;
	  }
	};
      };
      
      struct GlobalState
      {
        ALIGNED_CLASS;
      public:

        GlobalState () : numThreads(getNumberOfLogicalThreads()) {
	  threadStack = new WorkStack<BuildRecord,SIZE_WORK_STACK>[numThreads]; 
        }
        
        ~GlobalState () {
          delete[] threadStack;
        }

      public:
	size_t numThreads;
	WorkHeap<BuildRecord> heap;
        __aligned(64) WorkStack<BuildRecord,SIZE_WORK_STACK>* threadStack;
        ObjectPartition::ParallelBinner parallelBinner;
      };

    public:

      /*! Constructor. */
      BVH4BuilderFast (LockStepTaskScheduler* scheduler, BVH4* bvh, size_t listMode, size_t logBlockSize, size_t logSAHBlockSize, bool needVertices, size_t primBytes, const size_t minLeafSize, const size_t maxLeafSize);
      
      /*! Destructor */
      ~BVH4BuilderFast ();
   
      /* build function */
      virtual void build(size_t threadIndex, size_t threadCount);

      /* single threaded build */
      virtual void build_sequential(size_t threadIndex, size_t threadCount);

    public:
      TASK_SET_FUNCTION(BVH4BuilderFast,computePrimRefs);
      TASK_FUNCTION(BVH4BuilderFast,buildSubTrees);
      TASK_SET_FUNCTION(BVH4BuilderFast,build_parallel);

    public:

      /*! compute number of primitives */
      virtual size_t number_of_primitives() { return 0; } // FIXME: create base class without these functions
    
      /*! creates build primitive array (sequential version) */
      virtual void create_primitive_array_sequential(size_t threadIndex, size_t threadCount, PrimInfo& pinfo) {};
    
      /*! creates build primitive array (parallel version) */
      virtual void create_primitive_array_parallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimInfo& pinfo) {};
    
      /*! build mode */
      enum { RECURSE_SEQUENTIAL = 1, RECURSE_PARALLEL = 2, BUILD_TOP_LEVEL = 3 };
      
      /*! splitting function that selects between sequential and parallel mode */
      void split(BuildRecord& current, BuildRecord& left, BuildRecord& right, const size_t mode, const size_t threadID, const size_t numThreads);
      
      /*! perform sequential binning and splitting */
      void splitSequential(BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild, const size_t threadID, const size_t numThreads);
      
      /*! perform parallel binning and splitting */
      void splitParallel(BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild, const size_t threadID, const size_t threads);
      
      /*! creates a small leaf node */
      virtual void createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID) = 0;

      /*! creates a large leaf node */
      void createLeaf(BuildRecord& current, Allocator& nodeAlloc, Allocator& leafAlloc, size_t threadIndex, size_t threadCount);
      
      /*! select between recursion and stack operations */
      void recurse_continue(BuildRecord& current, Allocator& nodeAlloc, Allocator& leafAlloc, const size_t mode, const size_t threadID, const size_t numThreads);
      
      /*! recursive build function */
      void recurse(BuildRecord& current, Allocator& nodeAlloc, Allocator& leafAlloc, const size_t mode, const size_t threadID, const size_t numThreads);
      
      static void splitFallback(PrimRef * __restrict__ const primref, BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild);
    
    public:
      LockStepTaskScheduler* scheduler;
      std::unique_ptr<GlobalState> state;

      BVH4* bvh;                               //!< Output BVH
      size_t listMode;
      size_t logBlockSize;
      size_t logSAHBlockSize;
      size_t blocks(size_t N) { return (N+((1<<logBlockSize)-1)) >> logBlockSize; }
      bool needVertices;
      size_t primBytes; 
      size_t minLeafSize;
      size_t maxLeafSize;
      
    protected:
      TaskScheduler::Task task;
      
    public:
      PrimRef* prims;
      size_t bytesPrims;
      
    protected:
      size_t numPrimitives;
    };

    template<typename Primitive>
    class BVH4BuilderFastT : public BVH4BuilderFast
    {
    public:
      BVH4BuilderFastT (BVH4* bvh, Scene* scene, size_t listMode, size_t logBlockSize, size_t logSAHBlockSize, bool needVertices, size_t primBytes, const size_t minLeafSize, const size_t maxLeafSize, bool parallel);
      void createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID);
    public:
      Scene* scene;         //!< input scene
    };

    template<typename Primitive>
    class BVH4BezierBuilderFast : public BVH4BuilderFastT<Primitive>
    {
    public:
      BVH4BezierBuilderFast (BVH4* bvh, Scene* scene, size_t listMode);
      BVH4BezierBuilderFast (BVH4* bvh, BezierCurves* geom, size_t listMode);
      size_t number_of_primitives();
      void create_primitive_array_sequential(size_t threadIndex, size_t threadCount, PrimInfo& pinfo);
      void create_primitive_array_parallel  (size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimInfo& pinfo);
    public:
      BezierCurves* geom;   //!< input mesh
    };
    
    template<typename Primitive>
    class BVH4TriangleBuilderFast : public BVH4BuilderFastT<Primitive>
    {
    public:
      BVH4TriangleBuilderFast (BVH4* bvh, Scene* scene, size_t listMode);
      BVH4TriangleBuilderFast (BVH4* bvh, TriangleMesh* mesh, size_t listMode);
      size_t number_of_primitives();
      void create_primitive_array_sequential(size_t threadIndex, size_t threadCount, PrimInfo& pinfo);
      void create_primitive_array_parallel  (size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimInfo& pinfo) ;
    public:
      TriangleMesh* geom;   //!< input mesh
    };
  
    template<typename Primitive>
    class BVH4UserGeometryBuilderFastT : public BVH4BuilderFastT<Primitive>
    {
    public:
      BVH4UserGeometryBuilderFastT (BVH4* bvh, Scene* scene, size_t listMode);
      BVH4UserGeometryBuilderFastT (BVH4* bvh, UserGeometryBase* geom, size_t listMode);
      size_t number_of_primitives();
      void create_primitive_array_sequential(size_t threadIndex, size_t threadCount, PrimInfo& pinfo);
      void create_primitive_array_parallel  (size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimInfo& pinfo) ;
    public:
      UserGeometryBase* geom;   //!< input geometry
    };

    class BVH4TopLevelBuilderFastT : public BVH4BuilderFast
    {
    public:
      BVH4TopLevelBuilderFastT (LockStepTaskScheduler* scheduler, BVH4* bvh);
      void createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID);
      void build(size_t threadIndex, size_t threadCount, PrimRef* prims_i, size_t N)
      {
	this->prims_i = prims_i;
	this->N = N;
	BVH4BuilderFast::build(threadIndex,threadCount);
      }
      size_t number_of_primitives();
      void create_primitive_array_sequential(size_t threadIndex, size_t threadCount, PrimInfo& pinfo);
      void create_primitive_array_parallel  (size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimInfo& pinfo);
    public:
      PrimRef* prims_i;
      size_t N;
    };

    template<typename Primitive>
    class BVH4SubdivBuilderFast : public BVH4BuilderFastT<Primitive>
    {
    public:
      BVH4SubdivBuilderFast (BVH4* bvh, Scene* scene, size_t listMode);
      BVH4SubdivBuilderFast (BVH4* bvh, SubdivMesh* geom, size_t listMode);
      virtual void build(size_t threadIndex, size_t threadCount);

      size_t number_of_primitives();
      void create_primitive_array_sequential(size_t threadIndex, size_t threadCount, PrimInfo& pinfo);
      void create_primitive_array_parallel  (size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimInfo& pinfo);
    public:
      SubdivMesh* geom;   //!< input mesh
    };

    class BVH4SubdivQuadQuad4x4BuilderFast : public BVH4BuilderFastT<PrimRef>
    {
    public:
      BVH4SubdivQuadQuad4x4BuilderFast (BVH4* bvh, Scene* scene, size_t listMode);
      virtual void build(size_t threadIndex, size_t threadCount);

      size_t number_of_primitives();
      void create_primitive_array_sequential(size_t threadIndex, size_t threadCount, PrimInfo& pinfo);
      void create_primitive_array_parallel  (size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimInfo& pinfo);

      Scene::Iterator<SubdivMesh> iter;
      ParallelForForPrefixSumState<PrimInfo> pstate;
    };

    class BVH4SubdivGridBuilderFast : public BVH4BuilderFastT<PrimRef>
    {
    public:
      BVH4SubdivGridBuilderFast (BVH4* bvh, Scene* scene, size_t listMode);
      virtual void build(size_t threadIndex, size_t threadCount);

      size_t number_of_primitives();
      void create_primitive_array_sequential(size_t threadIndex, size_t threadCount, PrimInfo& pinfo);
      void create_primitive_array_parallel  (size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimInfo& pinfo);

      Scene::Iterator<SubdivMesh> iter;
      ParallelForForPrefixSumState<PrimInfo> pstate;
    };

    class BVH4SubdivGridEagerBuilderFast : public BVH4BuilderFastT<PrimRef>
    {
    public:
      BVH4SubdivGridEagerBuilderFast (BVH4* bvh, Scene* scene, size_t listMode);
      virtual void build(size_t threadIndex, size_t threadCount);

      size_t number_of_primitives();
      void create_primitive_array_sequential(size_t threadIndex, size_t threadCount, PrimInfo& pinfo);
      void create_primitive_array_parallel  (size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimInfo& pinfo);

      Scene::Iterator<SubdivMesh> iter;
      ParallelForForPrefixSumState<PrimInfo> pstate;
    };

    class BVH4SubdivGridLazyBuilderFast : public BVH4BuilderFastT<Grid::LazyLeaf*>
    {
    public:
      BVH4SubdivGridLazyBuilderFast (BVH4* bvh, Scene* scene, size_t listMode);
      virtual void build(size_t threadIndex, size_t threadCount);

      size_t number_of_primitives();
      void create_primitive_array_sequential(size_t threadIndex, size_t threadCount, PrimInfo& pinfo);
      void create_primitive_array_parallel  (size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimInfo& pinfo);

      Scene::Iterator<SubdivMesh> iter;
      ParallelForForPrefixSumState<PrimInfo> pstate;
    };

    class BVH4BuilderFastGeneric : public BVH4BuilderFast
    {
    public:
      struct MakeLeaf {
        virtual BVH4::NodeRef operator() (Allocator& alloc, const PrimRef* prims, size_t N) const = 0;
      };

    public:
      BVH4BuilderFastGeneric (BVH4* bvh, PrimRef* prims, size_t N, const MakeLeaf& makeLeaf, size_t listMode, 
                              size_t logBlockSize, size_t logSAHBlockSize, bool needVertices, size_t primBytes, const size_t minLeafSize, const size_t maxLeafSize);
      ~BVH4BuilderFastGeneric();
      virtual void build(size_t threadIndex, size_t threadCount);
      void createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID);

    public:
      size_t N;
      const MakeLeaf& makeLeaf;
    };


    class BVH4SubdivPatch1CachedBuilderFast : public BVH4BuilderFastT<PrimRef>
    {
    private:
      bool fastUpdateMode;
      size_t fastUpdateMode_numFaces;

      BBox3fa refit(NodeRef& ref);

    public:
      BVH4SubdivPatch1CachedBuilderFast (BVH4* bvh, Scene* scene, size_t listMode);
      virtual void build(size_t threadIndex, size_t threadCount);

      size_t number_of_primitives();
      void create_primitive_array_sequential(size_t threadIndex, size_t threadCount, PrimInfo& pinfo);
      void create_primitive_array_parallel  (size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimInfo& pinfo);

      void createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID);
      virtual void build_sequential(size_t threadIndex, size_t threadCount);

      Scene::Iterator<SubdivMesh> iter;
      ParallelForForPrefixSumState<PrimInfo> pstate;
    };

  }
}
