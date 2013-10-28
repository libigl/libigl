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

#ifndef __EMBREE_BVH4AOS_MIC_TASK_SCHEDULER_H__
#define __EMBREE_BVH4AOS_MIC_TASK_SCHEDULER_H__

#include "bvh4aos_globals.h"
#include "bvh4aos_builder_util.h"

namespace embree
{

  class QBVHTaskScheduler
  {
  public:

    static AtomicCounter taskCounter;
    static QuadTreeBarrier taskBarrier;

    static unsigned int NUM_TOTAL_THREADS;
    static unsigned int NUM_TOTAL_CORES;
    static unsigned int CONTROL_THREAD_ID;

    static AtomicMutex coreTaskCounter[MAX_MIC_CORES];
    static AtomicMutex debugMutex;

    static void (* taskPtr)(const unsigned int threadID);

    static bool taskDispatch(const unsigned int threadID = CONTROL_THREAD_ID);
    static void taskDispatchMainLoop(const unsigned int threadID);
    static void releaseThreads();
  
    static _INLINE bool dispatchTask(void (* task)(const unsigned int threadID),
				     const unsigned int threadID = CONTROL_THREAD_ID)
    {
      QBVHTaskScheduler::taskPtr = task;
      return QBVHTaskScheduler::taskDispatch(threadID);
    }

    static _INLINE void syncThreads(const unsigned int threadID,
				    const unsigned int threads = NUM_TOTAL_THREADS)
    {
      taskBarrier.sync(threadID,threads);
    }

  };

};

#endif
