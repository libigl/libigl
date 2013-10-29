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

#include "taskscheduler_mic.h"

namespace embree
{
  TaskSchedulerMIC::TaskSchedulerMIC() 
    : nextScheduleIndex(0)
  {
    for (size_t i=0; i<NUM_TASKS; i++) tasks[i] = NULL;
    for (size_t i=0; i<NUM_TASKS; i++) locks[i] = 0;
  }

  void TaskSchedulerMIC::add(ssize_t threadIndex, QUEUE queue, Task* task)
  {
    if (task->event) task->event->inc();
    size_t liveIndex = (nextScheduleIndex++)&(NUM_TASKS-1);
    if (tasks[liveIndex] || locks[liveIndex]) 
      throw std::runtime_error("task list full");
    __memory_barrier();
    locks[liveIndex] = numThreads;
    __memory_barrier();
    tasks[liveIndex] = task;
    __memory_barrier();
  }

  void TaskSchedulerMIC::run(size_t threadIndex, size_t threadCount)
  {
    size_t myIndex = 0;
    while (1) 
    {
      /* wait for available task */
      while (likely(!tasks[myIndex] && !terminateThreads)) {	  
        __pause(1023);
        continue;
      }
      
      /* terminate thread */
      if (terminateThreads) 
	{
	  return;
	}
	
      /* take next task from task list */
      Task* task = tasks[myIndex];
      TaskScheduler::Event* event = task->event;
      thread2event[threadIndex] = event; 

      while (true) 
      {
        ssize_t elt = --task->started;
        if (elt < 0) break;
        
        if (task->run) task->run(task->runData,threadIndex,threadCount,elt,task->elts,task->event);
      }
      
      /* free task slot */
      if (atomic_add(&locks[myIndex],-1) == 1) 
      {
        /* complete the task */
        if (task->complete) {
	  task->complete(task->completeData,threadIndex,threadCount,task->event);
	}
        if (event) event->dec();
        __memory_barrier();
        tasks[myIndex] = NULL;
      }

      /* goto next task slot */
      myIndex = (myIndex+1)&(NUM_TASKS-1);
    }
  }

  void TaskSchedulerMIC::terminate() {
    terminateThreads = true;
  }
}

