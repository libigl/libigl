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

#include "taskscheduler.h"
#include "sys/sync/mutex.h"
#include "sys/sync/barrier.h"
#include "sys/sync/atomic.h"

namespace embree
{
  /*! Task scheduler using active synchronization. */
  class TaskSchedulerMIC : public TaskScheduler
  {
  public:

    /* number of two, needs to be power of two */
    static const unsigned int NUM_TASKS = 1024;

    /*! construction */
    TaskSchedulerMIC();

  private:

    /*! adds a task to the specified task queue */
    void add(ssize_t threadIndex, QUEUE queue, Task* task);

    /*! thread function */
    void run(size_t threadIndex, size_t threadCount);

    /*! waits for an event out of a task */
    void wait(size_t threadIndex, size_t threadCount, Event* event);

    /*! processes next task */
    void work(size_t threadIndex, size_t threadCount,Task *task);

    /*! sets the terminate thread variable */
    void terminate();

  private:
    __aligned(64) AlignedAtomicCounter32 head_task_list; /*! next index in the task queue where we'll insert a live task */
    __aligned(64) Task* volatile tasks[NUM_TASKS]; //!< queue of tasks
    __aligned(64) LinearBarrierActive barrier;
  };
}
