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

#ifndef __EMBREE_TASKSCHEDULER_STANDARD_H__
#define __EMBREE_TASKSCHEDULER_STANDARD_H__

#include "taskscheduler.h"
#include "thread.h"
#include "sync/atomic.h"
#include "sync/mutex.h"
#include "sync/condition.h"

#include <vector>
#include <list>

namespace embree
{
  /*! Task scheduler implementing a stack of tasks. */
  class TaskSchedulerStandard : public TaskScheduler
  {
    /* initialization structure for threads */
    struct Thread {
      Thread (size_t tid, TaskSchedulerStandard* scheduler) : tid(tid), scheduler(scheduler) {}
      size_t tid;
      TaskSchedulerStandard* scheduler;
    };

    /* task class */
    class Task
    {
    public:
      __forceinline Task() 
        : run(NULL), runData(NULL), complete(NULL), completeData(NULL), name(NULL) {}

      __forceinline Task(runFunction run, void* runData, size_t elts, completeFunction complete, void* completeData, const char* name)
        : run(run), runData(runData), complete(complete), completeData(completeData), started(elts), completed(elts), name(name) {}
      
    public:
      runFunction run;             //!< run function
      void* runData;               //!< data pointer to execute run function
      completeFunction complete;   //!< complete function
      void* completeData;          //!< data pointer to execute complete function
      Atomic started;              //!< counts the number of started task set elements
      Atomic completed;            //!< counts the number of completed task set elements
      const char* name;            //!< name of this task
    };

  public:

    /*! task scheduler constructor */
    TaskSchedulerStandard (size_t numThreads = 0);

    /*! destructor */
    ~TaskSchedulerStandard();

    /*! returns number of threads created */
    size_t getNumThreads() const { return numThreads; }
    
    /*! before using the task system, one has to execute the start function */
    void start();

    /*! only a single initial task can be added from the main thread */
    void addTask(const ThreadInfo& thread, QUEUE queue,
                 runFunction run, void* runData, size_t elts = 1, 
                 completeFunction complete = NULL, void* completeData = NULL, 
                 const char* name = NULL);

    /*! the stop function waits until all tasks have finished executing */
    void stop();

  private:

    /*! thread function */
    void run(size_t tid);

    /*! thread function */
    static void threadFunction(Thread* thread);
    
    /* thread handling */
  private:
    volatile bool terminateThreads;
    std::vector<thread_t> threads;
    size_t numThreads;

  public:
    MutexSys mainMutex;           //!< the task system get locked while used

    /* stack of tasks */
  private:
    MutexSys taskMutex;           //!< protection of the task stack
    ConditionSys taskCondition;   //!< condition to signal if data is in task stack
    std::list<Task*> tasks;

    /* signaling finished */
  private:
    MutexSys finishMutex; 
    ConditionSys finishCondition;
    Atomic activeTasks;
  };
}

#endif

