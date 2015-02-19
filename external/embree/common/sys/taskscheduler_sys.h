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
#include "sys/sync/condition.h"

namespace embree
{
  /*! Task scheduler implementing a stack of tasks. */
  class __hidden TaskSchedulerSys : public TaskScheduler
  {
  public:

    /*! construction */
    TaskSchedulerSys();

  private:

    /*! only a single initial task can be added from the main thread */
    void add(ssize_t threadIndex, QUEUE queue, Task* task);

    /*! waits for an event out of a task */
    void wait(size_t threadIndex, size_t threadCount, Event* event);

    /*! processes next task */
    void work(size_t threadIndex, size_t threadCount, bool wait);

    /*! thread function */
    void run(size_t threadIndex, size_t threadCount);

    /*! sets the terminate thread variable */
    void terminate();
    
  private:
    MutexSys mutex;           //!< mutex to protect access to task list
    ConditionSys condition;   //!< condition to signal new tasks
    size_t begin,end;         //!< current range of tasks
    std::vector<Task*> tasks; //!< queue of tasks
  };
}
