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

#include "builders/primrefalloc.h"
#include "builders/primrefblock.h"
#include "geometry/primitive.h"
#include "builders/heuristic_split.h"
#include "builders/primrefgen.h"

namespace embree
{
  namespace isa
  {
    class BVH8Builder : public Builder
    {
      ALIGNED_CLASS;
    public:
      
      /*! Type shortcuts */
      typedef BVH8::Node    Node;
      typedef BVH8::NodeRef NodeRef;
      typedef atomic_set<PrimRefBlockT<PrimRef> > PrimRefList;
      typedef LinearAllocatorPerThread::ThreadAllocator Allocator;

      /*! the build record stores all information to continue the build of some subtree */
      struct BuildRecord 
      {
      public:
	__forceinline BuildRecord () {}

	__forceinline BuildRecord (size_t depth) : depth(depth), pinfo(empty) {}

	__forceinline BuildRecord (size_t depth, PrimRefList& prims, const PrimInfo& pinfo, const Split& split, NodeRef* dst)
	: depth(depth), prims(prims), pinfo(pinfo), split(split), dst(dst) {}
	
	__forceinline friend bool operator< (const BuildRecord& a, const BuildRecord& b) {
	  return a.pinfo.size() < b.pinfo.size();
	}
	
      public:
	NodeRef*   dst;      //!< Pointer to the parent node's reference to us
	size_t     depth;    //!< Depth of the root of this subtree.
	PrimRefList prims;    //!< The list of primitives.
	PrimInfo   pinfo;    //!< Bounding info of primitives.
	Split      split;    //!< The best split for the primitives.
      };
      
    public:
      
      /*! Constructor. */
      BVH8Builder (BVH8* bvh, Scene* scene, TriangleMesh* mesh, 
		    size_t mode, size_t logBlockSize, size_t logSAHBlockSize, float intCost, bool needVertices, 
		    size_t primBytes, const size_t minLeafSize, const size_t maxLeafSize);

      /*! Destructor*/
      ~BVH8Builder();

      /*! builder entry point */
      void build(size_t threadIndex, size_t threadCount);
   
      /*! build job */
      TASK_SET_FUNCTION(BVH8Builder,build_parallel);
      //void build_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount);
      
      /*! creates a leaf node */
      virtual NodeRef createLeaf(size_t threadIndex, Allocator& nodeAlloc, Allocator& leafAlloc, PrimRefList& prims, const PrimInfo& pinfo) = 0;
      
      /*! creates a large leaf by adding additional internal nodes */
      NodeRef createLargeLeaf(size_t threadIndex, Allocator& nodeAlloc, Allocator& leafAlloc, PrimRefList& prims, const PrimInfo& pinfo, size_t depth);
      
      /*! copies topmost nodes to improve memory layout */
      NodeRef layout_top_nodes(size_t threadIndex, Allocator& nodeAlloc, NodeRef node);

      /*! finds best possible split for a list of triangles */
      template<bool PARALLEL>
      const Split find(size_t threadIndex, size_t threadCount, size_t depth, PrimRefList& prims, const PrimInfo& pinfo, bool spatial);

      /*! creates a node from some build record */
      template<bool PARALLEL>
        static size_t createNode(size_t threadIndex, size_t threadCount, Allocator& nodeAlloc, Allocator& leafAlloc, 
                                 BVH8Builder* parent, BuildRecord& record, BuildRecord records_o[BVH8::N]);

      /*! continues build */
      void continue_build(size_t threadIndex, size_t threadCount, Allocator& nodeAlloc, Allocator& leafAlloc,
                          BuildRecord& record);

      /*! recursively finishes build */
      void finish_build(size_t threadIndex, size_t threadCount, Allocator& nodeAlloc, Allocator& leafAlloc,
                        BuildRecord& record);

      /*! performs some brute force restructuring of the tree */
      void restructureTree(NodeRef& ref, size_t depth);

    protected:
      Scene* scene;                       //!< input geometry
      TriangleMesh* mesh;                 //!< input triangle mesh
      PrimRefBlockAlloc<PrimRef> alloc;   //!< Allocator for primitive blocks
      BVH8* bvh;                          //!< Output BVH8
      LockStepTaskScheduler* scheduler;

      /*! build record task list */
    private:
      volatile atomic_t activeBuildRecords;
      MutexSys taskMutex;
      vector_t<BuildRecord> tasks;
     
      /*! builder configuration*/
    protected:
      size_t listMode;
      size_t minLeafSize;                 //!< minimal size of a leaf
      size_t maxLeafSize;                 //!< maximal size of a leaf
      bool enableSpatialSplits;
      size_t logSAHBlockSize;             //!< set to the logarithm of block size to use for SAH
      atomic_t remainingReplications;     //!< remaining replications allowed by spatial splits
      
      /*! primitive information */
    protected:
      size_t intCost;
      size_t logBlockSize;
      bool needVertices;
      size_t primBytes; 
      size_t blocks(size_t N) { return (N+((1<<logBlockSize)-1)) >> logBlockSize; }  
    };

    /*! specializes the builder for different leaf types */
    template<typename Triangle>
    class BVH8BuilderT : public BVH8Builder
    {
    public:
      BVH8BuilderT (BVH8* bvh, Scene* scene, size_t mode);
      BVH8BuilderT (BVH8* bvh, TriangleMesh* mesh, size_t mode);
      NodeRef createLeaf(size_t threadIndex, Allocator& nodeAlloc, Allocator& leafAlloc, PrimRefList& prims, const PrimInfo& pinfo);
    };
  }
}
