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

#include "sys/platform.h"
#include "sys/thread.h"
#include "sys/sync/event.h"
#include "sys/sync/atomic.h"
#include "sys/sync/barrier.h"
#include "sys/ref.h"

#include <vector>

namespace embree
{
  /*! Interface to different task scheduler implementations. */
  /* __hidden */
  class __hidden TaskScheduler : public RefCount
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
      AtomicCounter started;              //!< counts the number of started task set elements
      AtomicCounter completed;            //!< counts the number of completed task set elements
      const char* name;            //!< name of this task
      AtomicCounter locks;
    };

    /* an event that gets triggered by a task when completed */
    struct Event 
    {
      __forceinline Event (int activeTasks = 1, Event* other = NULL) 
        : activeTasks(activeTasks), other(other) 
      {
        if (other) other->inc();
      }

      void inc() { activeTasks++; }
      void dec() { 
        if (--activeTasks == 0) {
          Event* other = this->other;
          trigger(); // may cause this event to get destroyed
          if (other) other->dec();
        }
      }
      virtual void trigger() { --activeTasks; }
      bool triggered() const { return activeTasks == -1; }

    public:
      AtomicCounter activeTasks;  //!< number of tasks in flight
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

  protected:

    /*! construction */
    TaskScheduler();

  public:

    /*! single instance of task scheduler */
    static TaskScheduler* instance;
    
    /*! creates the threads */
    static void create(size_t numThreads = 0);

    /*! returns the number of threads used */
    static size_t getNumThreads();

    /*! enables specified number of threads */
    static size_t enableThreads(size_t N);

    /*! add a task to the scheduler */
    static void addTask(ssize_t threadIndex, QUEUE queue, Task* task);

    /*! executes a task, function returns if execution finished */
    static void executeTask(size_t threadIndex, size_t threadCount, runFunction run, void* runData, size_t elts, completeFunction complete, void* completeData, const char* name);

    static void executeTask(size_t threadIndex, size_t threadCount, runFunction run, void* runData, size_t elts, const char* name);

    static void executeTask(size_t threadIndex, size_t threadCount, completeFunction complete, void* completeData, const char* name);

    /*! destroys the task scheduler */
    static void destroy();

    /*! waits for an event out of a task */
    static void waitForEvent(Event* event); // use only from main thread !!!

      /*! enters lockstep taskscheduler, main thread returns false */
    static bool enter(size_t threadIndex, size_t threadCount);

    /*! leaves lockstep taskscheduler */
    static void leave(size_t threadIndex, size_t threadCount);
    
  protected:

    /*! creates all threads */
    void createThreads(size_t numThreads);

    /*! thread function */
    static void threadFunction(void* thread);

    /*! thread function */
    virtual void run(size_t threadIndex, size_t threadCount) = 0;

    /*! add a task */
    virtual void add(ssize_t threadIndex, QUEUE queue, Task* task) = 0;

    /*! waits for an event out of a task */
    virtual void wait(size_t threadIndex, size_t threadCount, Event* event) = 0;

    /*! sets the terminate thread variable */
    virtual void terminate() = 0;

    /*! destroys all threads */
    void destroyThreads();

    /* thread handling */
  protected:
    volatile bool terminateThreads;
    std::vector<thread_t> threads;
    bool defaultNumThreads;
    size_t numThreads;
    volatile size_t numEnabledThreads;
    __forceinline bool isEnabled(size_t threadIndex) const { return threadIndex < numEnabledThreads; }
  };

#define TASK_FUNCTION(Class,Name) \
  void Name (const size_t threadID, const size_t numThreads);           \
  static void task_##Name (void* data, const size_t threadID, const size_t numThreads) { \
    ((Class*)data)->Name(threadID,numThreads);                          \
  }

#define TASK_SET_FUNCTION(Class,Name) \
  void Name (const size_t threadID, const size_t numThreads, const size_t taskID, const size_t numTasks); \
  static void _##Name (void* data, const size_t threadID, const size_t numThreads, const size_t taskID, const size_t numTasks) { \
    ((Class*)data)->Name(threadID,numThreads,taskID,numTasks);			\
  }

#define LOCAL_TASK_FUNCTION(Class,Name) \
  void Name (const size_t localThreadID, const size_t globalThreadID);			\
  static void task_##Name (void* data, const size_t localThreadID, const size_t globalThreadID) { \
    ((Class*)data)->Name(localThreadID,globalThreadID);				\
  }
  

  class LockStepTaskScheduler
  {
  public:

    struct Init
    {
      Init (size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler)
      : threadIndex(threadIndex), threadCount(threadCount), scheduler(scheduler) 
      {
	if (threadCount) 
	{
	  if (threadIndex == 0) {
            scheduler->taskBarrier.init(threadCount);
            scheduler->threadCount = threadCount;
	    LockStepTaskScheduler::scheduler = scheduler;
          }
	  scheduler->barrier.wait(threadCount);
	  scheduler->enter(threadIndex,threadCount);
	}
      }
      
      ~Init () {
	if (threadIndex == 0 && threadCount != 0) {
	  scheduler->leave(threadIndex,threadCount);
	  scheduler->barrier.reset();
	  LockStepTaskScheduler::scheduler = NULL;
	}
      }

      size_t threadIndex, threadCount;
      LockStepTaskScheduler* scheduler;
    };

    static __thread LockStepTaskScheduler* scheduler;
    static LockStepTaskScheduler* instance();
    static void setInstance(LockStepTaskScheduler*);

    static const unsigned int CONTROL_THREAD_ID = 0;

    __aligned(64) BarrierActive barrier;

    typedef void (*runFunction)(void* data, const size_t threadID, const size_t numThreads);
    typedef void (*runFunction2)(void* data, const size_t threadID, const size_t numThreads, const size_t taskID, const size_t numTasks);

    __aligned(64) AlignedAtomicCounter32 taskCounter;
    __aligned(64) runFunction taskPtr;
    __aligned(64) runFunction2 taskPtr2;
    size_t numTasks;
    size_t numThreads; 

    size_t getNumThreads() const { return threadCount; }
    size_t threadCount;
    __aligned(64) void* volatile data;

#if defined(__MIC__)
    __aligned(64) QuadTreeBarrier taskBarrier;
    //__aligned(64) LinearBarrierActive taskBarrier; //FIXME

#else
    __aligned(64) Barrier taskBarrier;
#endif

    bool enter(size_t threadIndex, size_t threadCount);
    void leave(size_t threadIndex, size_t threadCount);

    bool dispatchTask(const size_t threadID, size_t numThreads);

    void dispatchTaskMainLoop(const size_t threadID, const size_t numThreads);
    void releaseThreads(const size_t numThreads);

    __forceinline bool dispatchTask(runFunction task, void* data) {
      return dispatchTask(task, data, 0, threadCount);
    }
  
    __forceinline bool dispatchTask(runFunction task, void* data, const size_t threadID, const size_t numThreads)
    {
      LockStepTaskScheduler::taskPtr = task;
      LockStepTaskScheduler::taskPtr2 = NULL;
      LockStepTaskScheduler::data = data;
      return LockStepTaskScheduler::dispatchTask(threadID, numThreads);
    }

    __forceinline bool dispatchTask(const size_t threadID, const size_t numThreads, runFunction2 task, void* data, const size_t numTasks, const char* name = "")
    {
      LockStepTaskScheduler::taskPtr = NULL;
      LockStepTaskScheduler::taskPtr2 = task;
      LockStepTaskScheduler::data = data;
      LockStepTaskScheduler::numTasks = numTasks;
      return LockStepTaskScheduler::dispatchTask(threadID, numThreads);
    }

    __forceinline bool dispatchTaskSet(runFunction2 task, void* data, const size_t numTasks)
    {
      LockStepTaskScheduler::taskPtr = NULL;
      LockStepTaskScheduler::taskPtr2 = task;
      LockStepTaskScheduler::data = data;
      LockStepTaskScheduler::numTasks = numTasks;
      return LockStepTaskScheduler::dispatchTask(0, threadCount);
    }

    void syncThreads(const size_t threadID, const size_t numThreads);

    void syncThreadsWithReduction(const size_t threadID, 
				  const size_t numThreads,
				  void (* reductionFct)(const size_t currentThreadID,
							const size_t childThreadID,
							void *ptr),
				  void *ptr);



  };

  extern __aligned(64) LockStepTaskScheduler g_regression_task_scheduler;
  
  class __aligned(64) LockStepTaskScheduler4ThreadsLocalCore
  {
  public:

    static const unsigned int WAIT_CYCLES       = 256;
    
    void (* taskPtr)(void* data, const size_t localThreadID, const size_t globalThreadID);
    void* volatile data;
    volatile unsigned char threadState[2][4];
    volatile unsigned int mode;

    __aligned(64) AtomicMutex mutex;


    LockStepTaskScheduler4ThreadsLocalCore();
    bool dispatchTask(const size_t localThreadID, const size_t globalThreadID);
    
    void dispatchTaskMainLoop(const size_t localThreadID, const size_t globalThreadID);
    void releaseThreads(const size_t localThreadID, const size_t globalThreadID);
    
    __forceinline bool dispatchTask(void (* task)(void* data, const size_t localThreadID, const size_t globalThreadID),
				    void* _data, 
				    const size_t localThreadID,
				    const size_t globalThreadID)
    {
      taskPtr = task;
      data    = _data;
      return dispatchTask(localThreadID,globalThreadID);
    }
    
    void syncThreads(const size_t localThreadID);
  };
}
