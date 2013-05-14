// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#ifndef __EMBREE_TASKSCHEDULER_H__
#define __EMBREE_TASKSCHEDULER_H__

#include "sys/platform.h"
#include "ref.h"

namespace embree
{
  /*! Interface to different task scheduler implementations. */
  class TaskScheduler : public RefCount
  {
  public:

    /*! Task queues */
    enum QUEUE { GLOBAL_FRONT, GLOBAL_BACK };

    /*! thread and core ID struct */
    struct ThreadInfo 
    {
      ThreadInfo () : id(0) {}
      ThreadInfo (size_t id) : id(id) {}
      size_t id;  // thread ID from 0 to getNumberOfLogicalThreads()-1
    };
    
    /*! the run function executed for each work item of the task */
    typedef void (*runFunction)(const ThreadInfo& thread, void* ptr, size_t elt);
    
    /*! complete function executed at the end of the task */
    typedef void (*completeFunction)(const ThreadInfo& thread, void* ptr);

    /*! returns the number of threads used */
    virtual size_t getNumThreads() const = 0;

    /*! starts the tasking system */
    virtual void start() = 0;

    /*! stops the tasking system and waits for all tasks to finish */
    virtual void stop() = 0;

    /*! adds a task */
    virtual void addTask(const ThreadInfo& thread, QUEUE queue,
                         runFunction run, void* runData, size_t elts = 1, 
                         completeFunction complete = NULL, void* completeData = NULL, 
                         const char* name = NULL) = 0;
  };

  extern TaskScheduler* scheduler;
}

#endif

