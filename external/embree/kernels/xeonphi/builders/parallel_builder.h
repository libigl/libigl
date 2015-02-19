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

#include "common/alloc.h"
#include "common/accel.h"
#include "common/scene.h"
#include "geometry/primitive.h"
#include "builders/builder_util.h"
#include "builders/binning.h"

namespace embree
{
  class ParallelBuilderInterface : public Builder
  {
    ALIGNED_CLASS;
  public:
    static const size_t ALLOCATOR_NODE_BLOCK_SIZE = 64;
    typedef AtomicIDBlock<ALLOCATOR_NODE_BLOCK_SIZE> NodeAllocator;    

    /*! build mode */
    enum { RECURSE = 1, FILL_LOCAL_QUEUES = 2, BUILD_TOP_LEVEL = 3 };

    ParallelBuilderInterface(void* geometry)
    : scene((Scene*)geometry), 
      enablePerCoreWorkQueueFill(true),
      enableTaskStealing(true),
      numPrimitives((size_t)-1),
      atomicID(0),
      numNodes(0),
      numAllocated64BytesBlocks(0)
	{
    
	}

    virtual void build            (const size_t threadIndex, const size_t threadCount) = 0;
    virtual void allocateData     (const size_t threadCount, const size_t newNumPrimitives) = 0;
    virtual void computePrimRefs  (const size_t threadIndex, const size_t threadCount) = 0;
    virtual void createAccel      (const size_t threadIndex, const size_t threadCount) = 0;

    virtual size_t getNumPrimitives() = 0;
    virtual void printBuilderName()   = 0;

    virtual void buildSubTree(BuildRecord& current, 
			      NodeAllocator& alloc, 
			      const size_t mode,
			      const size_t threadID, 
			      const size_t numThreads) = 0;


    virtual std::string getStatistics() = 0;

  protected:
    Scene* scene;                 //!< input geometry
    size_t numPrimitives;
    size_t numNodes;
    size_t numAllocated64BytesBlocks;

    /*! bounds shared among threads */    
    __aligned(64) Centroid_Scene_AABB global_bounds;

    /*! global node allocator */
    __aligned(64) AlignedAtomicCounter32  atomicID;

    static const size_t SIZE_GLOBAL_WORK_STACK = 512;
    static const size_t SIZE_LOCAL_WORK_STACK  = 16;

    /*! global work queue */
    __aligned(64) WorkStack<BuildRecord,SIZE_GLOBAL_WORK_STACK> global_workStack;

    /*! local per core work queue */
    __aligned(64) WorkStack<BuildRecord,SIZE_LOCAL_WORK_STACK> local_workStack[MAX_MIC_CORES];

    /*! per core lock-step task scheduler */
    __aligned(64) LockStepTaskScheduler4ThreadsLocalCore localTaskScheduler[MAX_MIC_CORES];

    /*! flags to enable/disable per core work queues and task stealing */
    bool enablePerCoreWorkQueueFill;
    bool enableTaskStealing;

    /*! parallel task functions */
    TASK_FUNCTION(ParallelBuilderInterface,fillLocalWorkQueues);
    TASK_FUNCTION(ParallelBuilderInterface,buildSubTrees);    
  };


 class ParallelBinnedSAHBuilder : public ParallelBuilderInterface
  {
    ALIGNED_CLASS;
  public:

    ParallelBinnedSAHBuilder (void* geometry) : ParallelBuilderInterface(geometry) {}

  protected:

    /* single structure for all worker threads */
    __aligned(64) SharedBinningPartitionData global_sharedData;

    /* one 16-bins structure per thread */
    __aligned(64) Bin16 global_bin16[MAX_MIC_THREADS];

    /* one shared binning/partitoning structure per core */
    __aligned(64) SharedBinningPartitionData local_sharedData[MAX_MIC_CORES];


    static void reduceBinsParallel(const size_t currentThreadID,
				   const size_t childThreadID,
				   void *ptr);

  };

};
