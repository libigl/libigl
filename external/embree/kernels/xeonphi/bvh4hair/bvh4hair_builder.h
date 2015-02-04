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

#include "builders/parallel_builder.h"
#include "builders/builder_util.h"
#include "builders/binning.h"

#include "bvh4i/bvh4i_builder.h"
#include "bvh4hair.h"
#include "geometry/bezier1i.h"


namespace embree
{

  /*! derived binned-SAH builder supporting hair primitives */  
  class BVH4HairBuilder : public ParallelBinnedSAHBuilder
  {
    ALIGNED_CLASS;
    
  protected:
    static const size_t ALLOCATOR_NODE_BLOCK_SIZE = 64;
    typedef AtomicIDBlock<ALLOCATOR_NODE_BLOCK_SIZE> NodeAllocator;    

    static const size_t MAX_ITEMS_PER_LEAF = 2;

  public:

    class __aligned(64) BuildRecordOBB : public BuildRecord
    {
    public:
      BuildRecordOBB() {}

      BuildRecordOBB(const BuildRecord &b) 
	{
	  *(BuildRecord*)this = b;
	  xfm = LinearSpace3fa( zero );
	}

      __aligned(64) LinearSpace3fa xfm;

      __forceinline void PreQuantizeMatrix()
      {
	mic_f col0 = broadcast4to16f(&xfm.vx) * 127.0f;
	mic_f col1 = broadcast4to16f(&xfm.vy) * 127.0f;
	mic_f col2 = broadcast4to16f(&xfm.vz) * 127.0f;
	
	mic_i char_col0,char_col1,char_col2;
	store16f_int8(&char_col0,col0);
	store16f_int8(&char_col1,col1);
	store16f_int8(&char_col2,col2);

	mic_f new_col0 = load16f_int8((char*)&char_col0);
	mic_f new_col1 = load16f_int8((char*)&char_col1);
	mic_f new_col2 = load16f_int8((char*)&char_col2);

	store4f(&xfm.vx,new_col0);
	store4f(&xfm.vy,new_col1);
	store4f(&xfm.vz,new_col2);

      }

      __forceinline friend std::ostream &operator<<(std::ostream &o, const BuildRecordOBB &br)
	{
	  o << "centroid2 = " << br.bounds.centroid2 << " ";
	  o << "geometry  = " << br.bounds.geometry << " ";
	  o << "begin       " << br.begin << " ";
	  o << "end         " << br.end << " ";
	  o << "items       " << br.end-br.begin << " ";
	  o << "parentPtr   " << br.parentPtr << " ";
	  o << "flags       " << br.flags << " ";
	  o << "sArea       " << br.sArea << " ";
	  o << "matrix      " << br.xfm << " ";
	  return o;
	};

    };

    BVH4Hair *bvh4hair;
    Bezier1i *prims;
    Bezier1i *accel;
    BVH4Hair::UnalignedNode*   node;
    
    size_t size_prims;
    size_t size_accel;
    size_t size_nodes;

    const size_t num64BytesBlocksPerNode;
    
  BVH4HairBuilder(BVH4Hair* bvh, void* geometry) 
    : ParallelBinnedSAHBuilder(geometry),
      bvh4hair(bvh),
      prims(NULL),
      node(NULL),
      accel(NULL),
      size_prims(0),
      size_nodes(0),
      size_accel(0),
      num64BytesBlocksPerNode(4)
	{
	}

    virtual ~BVH4HairBuilder() 
      {	
      }

    virtual size_t getNumPrimitives();
    virtual void   printBuilderName();
    virtual void   allocateData   (const size_t threadCount, const size_t newNumPrimitives);
    virtual void   computePrimRefs(const size_t threadIndex, const size_t threadCount);
    virtual void   build          (const size_t threadIndex, const size_t threadCount);
    virtual void   createAccel    (const size_t threadIndex, const size_t threadCount);

    virtual void buildSubTree(BuildRecord& current, 
			      NodeAllocator& alloc, 
			      const size_t mode,
			      const size_t threadID, 
			      const size_t numThreads);

    virtual std::string getStatistics();

    void build_main(size_t threadIndex, size_t threadCount);

  protected:

    __forceinline void createBVH4HairLeaf(void *parentPtr,
					  unsigned int offset, 
					  unsigned int items)
    {
      assert(items <= BVH4Hair::N);
      const unsigned int v = (offset << BVH4Hair::encodingBits) | BVH4Hair::leaf_mask | (items-1);
      unsigned int *ptr = (unsigned int*)parentPtr;
      *ptr = v;
    }

    __forceinline void createBVH4HairNode(void *parentPtr,
					  unsigned int index,
					  const size_t flags = 0)
    {
      unsigned int *ptr = (unsigned int*)parentPtr;
      *ptr = (index*sizeof(BVH4Hair::UnalignedNode)) | flags;      
    }


    /*! splitting function that selects between sequential and parallel mode */
    bool split(BuildRecord& current, BuildRecord& left, BuildRecord& right, const size_t mode, const size_t threadID, const size_t numThreads);

    /*! perform sequential binning and splitting */
    bool splitSequential(BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild);

    /*! perform parallel binning and splitting using all threads on all cores*/
    bool splitParallelGlobal(BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild, const size_t threadID, const size_t threads);

    /*! perform parallel binning and splitting using only the threads per core*/
    bool splitParallelLocal(BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild, const size_t threadID);

    /*! perform parallel splitting */
    void parallelPartitioning(BuildRecord& current,
			      Bezier1i * __restrict__ l_source,
			      Bezier1i * __restrict__ r_source,
			      Bezier1i * __restrict__ l_dest,
			      Bezier1i * __restrict__ r_dest,
			      const Split &split,
			      Centroid_Scene_AABB &local_left,
			      Centroid_Scene_AABB &local_right);			      

    void allocateMemoryPools(const size_t numPrims, const size_t numNodes);

    /*! recursive build functions */
    void recurseSAH(BuildRecord& current, NodeAllocator& alloc, const size_t mode, const size_t threadID, const size_t numThreads);
    void createLeaf(BuildRecord& current, NodeAllocator& alloc,const size_t threadIndex, const size_t threadCount);
    void recurse(BuildRecord& current, NodeAllocator& alloc,const size_t mode, const size_t threadID, const size_t numThreads);

    /*! unaligned splits */
    void recurseOBB(BuildRecordOBB& current, NodeAllocator& alloc, const size_t mode, const size_t threadID, const size_t numThreads);

    /*! perform sequential binning and splitting */
    bool splitSequentialOBB(BuildRecordOBB& current, BuildRecordOBB& leftChild, BuildRecordOBB& rightChild, const bool binAABB = false);

    void computeUnalignedSpace( BuildRecordOBB& current );
    void computeUnalignedSpaceBounds( BuildRecordOBB& current );

    TASK_FUNCTION(BVH4HairBuilder,computePrimRefsBezierCurves);
    TASK_FUNCTION(BVH4HairBuilder,parallelBinningGlobal);
    TASK_FUNCTION(BVH4HairBuilder,parallelPartitioningGlobal);
    LOCAL_TASK_FUNCTION(BVH4HairBuilder,parallelBinningLocal);
    LOCAL_TASK_FUNCTION(BVH4HairBuilder,parallelPartitioningLocal);

  };



}
