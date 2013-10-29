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

#ifndef __EMBREE_TASKSCHEDULER_MIC_H__
#define __EMBREE_TASKSCHEDULER_MIC_H__

#include "taskscheduler.h"
#include "sys/sync/mutex.h"

namespace embree
{
  /*! Task scheduler using active synchronization. */
  class TaskSchedulerMIC : public TaskScheduler
  {
  public:

    enum { NUM_TASKS = 4*1024 };

    /*! construction */
    TaskSchedulerMIC();

  private:

    /*! adds a task to the specified task queue */
    void add(ssize_t threadIndex, QUEUE queue, Task* task);

    /*! thread function */
    void run(size_t threadIndex, size_t threadCount);

    /*! sets the terminate thread variable */
    void terminate();

  private:
    Atomic nextScheduleIndex; /*! next index in the task queue where we'll insert a live task */
    Task* volatile tasks[NUM_TASKS]; //!< queue of tasks
    volatile atomic_t locks[NUM_TASKS]; 
  };
}

#endif

