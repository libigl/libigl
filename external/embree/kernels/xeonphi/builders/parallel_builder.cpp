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

#include "parallel_builder.h"

namespace embree
{

  static double dt = 0.0f;

  void ParallelBuilderInterface::fillLocalWorkQueues(const size_t threadID, const size_t numThreads)
  {
    __aligned(64) BuildRecord br;
    const size_t numCores = (numThreads+3)/4;
    const size_t coreID   = threadID/4;

    if (threadID % 4 == 0)
      {
	unsigned int ID = coreID; // LockStepTaskScheduler::taskCounter.inc();
	
	while (true) 
	  {
	    /* get build record from global queue */
	    if (ID >= global_workStack.size()) break;
	    br = global_workStack.get(ID);
	    bool success = local_workStack[coreID].push(br);	  
	    if (!success) FATAL("can't fill local work queues");
	    ID += numCores;
	  }    

      }
  }

  void ParallelBuilderInterface::buildSubTrees(const size_t threadID, const size_t numThreads)
  {
    NodeAllocator alloc(atomicID,numAllocated64BytesBlocks);
    __aligned(64) BuildRecord br;
    const size_t numCores = (numThreads+3)/4;
    const size_t globalCoreID   = threadID/4;

    const size_t MIN_LOCAL_WORK_QUEUE_ENTRIES = 8;

    if (enablePerCoreWorkQueueFill && numThreads > 1)
      {
	const size_t globalThreadID = threadID;
	const size_t localThreadID  = threadID % 4;
    
	if (localThreadID != 0)
	  {
	    localTaskScheduler[globalCoreID].dispatchTaskMainLoop(localThreadID,globalThreadID);
	  }
	else
	  {
	    local_workStack[globalCoreID].mutex.inc();
	    while (local_workStack[globalCoreID].size() < MIN_LOCAL_WORK_QUEUE_ENTRIES && 
		   local_workStack[globalCoreID].size()+ 4 <= SIZE_LOCAL_WORK_STACK) 
	      {
		BuildRecord br;
		if (!local_workStack[globalCoreID].pop_largest(br)) break;
		buildSubTree(br,alloc,FILL_LOCAL_QUEUES,globalThreadID,4);
	      }

	    localTaskScheduler[globalCoreID].releaseThreads(localThreadID,globalThreadID);	
	    local_workStack[globalCoreID].mutex.dec();
	  }
      }

    while(true)
      {
      /* process local work queue */
	while (1)
	  {
	    if (!local_workStack[globalCoreID].pop_largest(br)) 
	      {
		if (local_workStack[globalCoreID].mutex.val() > 0)
		  {
		    __pause_cpu(1024);
		    continue;
		  }
		else
		  break;
	      }
	    local_workStack[globalCoreID].mutex.inc();
	    buildSubTree(br,alloc,RECURSE,threadID,numThreads);
	    local_workStack[globalCoreID].mutex.dec();
	  }

	/* try task stealing */
        bool success = false;
	if (enableTaskStealing && numThreads > 4)
	  {
	    for (size_t i=0; i<numThreads; i++)
	      {
		unsigned int next_threadID = (threadID+i);
		if (next_threadID >= numThreads) next_threadID -= numThreads;
		const unsigned int next_globalCoreID   = next_threadID/4;

		assert(next_globalCoreID < numCores);
		if (local_workStack[next_globalCoreID].pop_smallest(br)) { 
		  success = true;
		  break;
		}
	      }
	  }
        if (!success) break; 
	
	local_workStack[globalCoreID].mutex.inc();
	buildSubTree(br,alloc,RECURSE,threadID,numThreads);
	local_workStack[globalCoreID].mutex.dec();

      }

  }


  void ParallelBinnedSAHBuilder::reduceBinsParallel(const size_t currentThreadID,
						    const size_t childThreadID,
						    void *ptr)
  {
    Bin16 *__restrict__ bin16 = (Bin16*)ptr;
    bin16[childThreadID].prefetchL2();
    bin16[currentThreadID].merge(bin16[childThreadID]);
  }



};
