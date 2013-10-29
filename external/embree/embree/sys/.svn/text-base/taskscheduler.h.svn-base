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

#ifndef __EMBREE_TASKSCHEDULER_H__
#define __EMBREE_TASKSCHEDULER_H__

#include "sys/platform.h"
#include "sys/thread.h"
#include "sys/sync/event.h"
#include "sys/sync/atomic.h"
#include "sys/ref.h"

#include <vector>

namespace embree
{
  /*! Interface to different task scheduler implementations. */
  class TaskScheduler : public RefCount
  {
  public:
    struct Event;

    /*! Task queues */
    enum QUEUE { GLOBAL_FRONT, GLOBAL_BACK };

#define TASK_RUN_FUNCTION(Class,name)                                   \
    void name(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* taskGroup); \
    static void _##name(void* This, size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* taskGroup) { \
      ((Class*)This)->name(threadIndex,threadCount,taskIndex,taskCount,taskGroup); \
    }
    
#define TASK_COMPLETE_FUNCTION(Class,name)                              \
    void name(size_t threadIndex, size_t threadCount, TaskScheduler::Event* taskGroup); \
    static void _##name(void* This, size_t threadIndex, size_t threadCount, TaskScheduler::Event* taskGroup) { \
      ((Class*)This)->name(threadIndex,threadCount,taskGroup); \
    }

#define TASK_RUN_FUNCTION_(Class,name)                                   \
    static void _##name(void* This, size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* taskGroup) { \
      ((Class*)This)->name(threadIndex,threadCount,taskIndex,taskCount,taskGroup); \
    }
    
#define TASK_COMPLETE_FUNCTION_(Class,name)                              \
    static void _##name(void* This, size_t threadIndex, size_t threadCount, TaskScheduler::Event* taskGroup) { \
      ((Class*)This)->name(threadIndex,threadCount,taskGroup); \
    }
    
    /*! the run function executed for each work item of the task */
    typedef void (*runFunction)(void* data, size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, Event* taskGroup);
    
    /*! complete function executed at the end of the task */
    typedef void (*completeFunction)(void* data, size_t threadIndex, size_t threadCount, Event* taskGroup);

    /* task class */
    class Task
    {
    public:
      __forceinline Task() 
        : event(NULL), run(NULL), runData(NULL), complete(NULL), completeData(NULL), name(NULL), locks(0) {}

      __forceinline Task(Event* event, runFunction run, void* runData, size_t elts, completeFunction complete, void* completeData, const char* name)
        : event(event), run(run), runData(runData), elts(elts), complete(complete), completeData(completeData), 
        started(elts), completed(elts), name(name), locks(0) {}

      __forceinline Task(Event* event, completeFunction complete, void* completeData, const char* name)
        : event(event), run(NULL), runData(NULL), elts(1), complete(complete), completeData(completeData), 
        started(1), completed(1), name(name), locks(0) {}

    public:
      Event* event;
      runFunction run;             //!< run function
      void* runData;               //!< data pointer to execute run function
      size_t elts;                 //!< total number of elements
      completeFunction complete;   //!< complete function
      void* completeData;          //!< data pointer to execute complete function
      Atomic started;              //!< counts the number of started task set elements
      Atomic completed;            //!< counts the number of completed task set elements
      const char* name;            //!< name of this task
      Atomic locks;
    };

    /* an event that gets triggered by a task when completed */
    struct Event 
    {
      __forceinline Event() {}
      __forceinline Event (int activeTasks, Event* other) 
        : activeTasks(activeTasks), other(other) 
      {
        if (other) other->inc();
      }

      __forceinline void inc() { activeTasks++; }
      __forceinline void dec() { 
        if (--activeTasks == 0) {
          Event* other = this->other;
          trigger(); // may cause this event to get destroyed
          if (other) other->dec();
        }
      }
      virtual void trigger() = 0;

    public:
      Atomic activeTasks;  //!< number of tasks in flight
      Event* other;
    };

    /* a group of tasks that one can wait for */
    struct EventSync : public Event
    {
      __forceinline EventSync (Event* other = NULL) : Event(1,other) {}
      __forceinline void sync() { dec(); event.wait(); }
      void trigger() { event.signal(); }
    public:
      EventSys event;
    };

    /* group of tasks that calls a completion function */
    struct EventScheduleTask : public Event
    {
      __forceinline EventScheduleTask() {}
      __forceinline EventScheduleTask (Event* other, completeFunction complete, void* completeData, const char* name) 
        : Event(0,other), task(other,NULL,NULL,1,complete,completeData,name) {}
      void trigger() { 
        addTask(-1,GLOBAL_FRONT,&task);
      }
    public:
      Task task;
    };

  protected:

    /*! construction */
    TaskScheduler();

  public:

    /*! single instance of task scheduler */
    static TaskScheduler* instance;

    /*! creates the threads */
    static void create(size_t numThreads);

    /*! returns the number of threads used */
    static size_t getNumThreads();

    /*! add a task to the scheduler */
    static void addTask(ssize_t threadIndex, QUEUE queue, Task* task);

    /*! destroys the task scheduler */
    static void destroy();

    /*! returns ISPC event of the thread */
    static Event* getISPCEvent(ssize_t threadIndex);
    
  protected:

    /*! creates all threads */
    void createThreads(size_t numThreads);

    /*! thread function */
    static void threadFunction(void* thread);

    /*! thread function */
    virtual void run(size_t threadIndex, size_t threadCount) = 0;

    /*! add a task */
    virtual void add(ssize_t threadIndex, QUEUE queue, Task* task) = 0;

    /*! sets the terminate thread variable */
    virtual void terminate() = 0;

    /*! destroys all threads */
    void destroyThreads();

    /* thread handling */
  protected:
    volatile bool terminateThreads;
    std::vector<thread_t> threads;
    size_t numThreads;
    Event** thread2event;
  };
}

#endif

