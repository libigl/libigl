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

#include "bvh4aos_task_scheduler.h"

#define DBG_THREADS(x) 
#define OUTPUT(x) 

namespace embree
{
  using namespace std;

  MIC_ALIGN AtomicCounter QBVHTaskScheduler::taskCounter           = 0;
  unsigned int  QBVHTaskScheduler::NUM_TOTAL_THREADS               = 1;
  unsigned int  QBVHTaskScheduler::NUM_TOTAL_CORES                 = 1;
  unsigned int  QBVHTaskScheduler::CONTROL_THREAD_ID               = 0;

  MIC_ALIGN void (* QBVHTaskScheduler::taskPtr)(const unsigned int threadID) = NULL;
  MIC_ALIGN QuadTreeBarrier QBVHTaskScheduler::taskBarrier;
  MIC_ALIGN AtomicMutex QBVHTaskScheduler::debugMutex;

  void QBVHTaskScheduler::taskDispatchMainLoop(const unsigned int threadID)
  {
    while(1)	 
      {
	bool dispatch = taskDispatch(threadID);
	if (dispatch == true) break; 
      }  
  }

  bool QBVHTaskScheduler::taskDispatch(const unsigned int threadID)
  {
    // =========================================================
    if (threadID == CONTROL_THREAD_ID)
      {
	DBG_THREADS(cout << "reset task counter " << endl << flush);
	QBVHTaskScheduler::taskCounter.reset(0);
      }

    DBG_THREADS(cout << "enter start barrier " << threadID << endl << flush);
    syncThreads(threadID,NUM_TOTAL_THREADS);
    DBG_THREADS(cout << "leave start barrier " << threadID << endl << flush);
  

    // =========================================================

    if (taskPtr == NULL) 
      {
	return true;
      }
    else
      {
	DBG_THREADS(cout << "executing task code for thread " << threadID << endl << flush);
	(*taskPtr)(threadID); // === jump to task code ===
      }

    // =========================================================
    DBG_THREADS(cout << "enter end barrier " << threadID << endl << flush);
    syncThreads(threadID,NUM_TOTAL_THREADS);
    DBG_THREADS(cout << "leave end barrier " << threadID << endl << flush);


    // =========================================================

    return false;
  }


  void QBVHTaskScheduler::releaseThreads()
  {
    taskPtr = NULL;
    QBVHTaskScheduler::taskDispatch(CONTROL_THREAD_ID);  
  }

};
