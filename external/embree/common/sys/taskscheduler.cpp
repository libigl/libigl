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

#include "taskscheduler.h"
#include "taskscheduler_sys.h"
#if defined(__MIC__)
#include "taskscheduler_mic.h"
#endif
#include "sysinfo.h"
#include "tasklogger.h"
#include "sys/sync/atomic.h"

#define DBG_THREADS(x)

namespace embree
{
  /* initialization structure for threads */
  struct Thread 
  {
    Thread (size_t threadIndex, size_t threadCount, TaskScheduler* scheduler) 
      : threadIndex(threadIndex), threadCount(threadCount), scheduler(scheduler) {}

  public:
    size_t threadIndex;
    size_t threadCount;
    TaskScheduler* scheduler;
  };
  
  TaskScheduler* TaskScheduler::instance = NULL;

  void TaskScheduler::create(size_t numThreads)
  {
    if (instance)
      throw std::runtime_error("Embree threads already running.");

    /* enable fast pthreads tasking system */
#if defined(__MIC__)
    instance = new TaskSchedulerMIC; 
    //instance = new TaskSchedulerSys; 
#else
    instance = new TaskSchedulerSys; 
#endif

#if 1
    instance->createThreads(numThreads);
#else
    instance->createThreads(1);
    std::cout << "WARNING: Using only a single thread." << std::endl;
#endif
  }

  size_t TaskScheduler::getNumThreads() 
  {
    if (!instance) throw std::runtime_error("Embree tasks not running.");
    return instance->numThreads;
  }

  void TaskScheduler::addTask(ssize_t threadIndex, QUEUE queue, Task* task)
  {
    if (!instance) throw std::runtime_error("Embree tasks not running.");
    instance->add(threadIndex,queue,task);
  }

  void TaskScheduler::executeTask(size_t threadIndex, size_t threadCount, 
                                  runFunction run, void* runData, size_t elts, completeFunction complete, void* completeData, const char* name)
  {
    TaskScheduler::Event event;
    TaskScheduler::Task task(&event,run,runData,elts,complete,completeData,name);
    instance->add(threadIndex,TaskScheduler::GLOBAL_FRONT,&task);
    instance->wait(threadIndex,threadCount,&event);
  }

  void TaskScheduler::executeTask(size_t threadIndex, size_t threadCount,  
                                  runFunction run, void* runData, size_t elts, const char* name)
  {
    TaskScheduler::Event event;
    TaskScheduler::Task task(&event,run,runData,elts,NULL,NULL,name);
    instance->add(threadIndex,TaskScheduler::GLOBAL_FRONT,&task);
    instance->wait(threadIndex,threadCount,&event);
  }

  void TaskScheduler::executeTask(size_t threadIndex, size_t threadCount,  
                                  completeFunction complete, void* completeData, const char* name)
  {
    TaskScheduler::Event event;
    TaskScheduler::Task task(&event,NULL,NULL,1,complete,completeData,name);
    instance->add(threadIndex,TaskScheduler::GLOBAL_FRONT,&task);
    instance->wait(threadIndex,threadCount,&event);
  }


  void TaskScheduler::waitForEvent(Event* event) {
    instance->wait(0,instance->getNumThreads(),event);
  }

  void TaskScheduler::destroy() 
  {
    if (instance) {
      instance->destroyThreads();
      delete instance; 
      instance = NULL;
    }
  }
  
  TaskScheduler::Event* TaskScheduler::getISPCEvent(ssize_t threadIndex)
  {
    if (!instance) throw std::runtime_error("Embree tasks not running.");
    if (threadIndex < 0 || threadIndex >= (ssize_t)instance->numThreads)
      throw std::runtime_error("invalid thread index");

    return instance->thread2event[threadIndex].event;
  }

  TaskScheduler::TaskScheduler () 
    : terminateThreads(false), numThreads(0), thread2event(NULL) {}

  void TaskScheduler::createThreads(size_t numThreads_in)
  {
    numThreads = numThreads_in;
#if defined(__MIC__)
    if (numThreads == 0) numThreads = getNumberOfLogicalThreads()-4;
#else
    if (numThreads == 0) numThreads = getNumberOfLogicalThreads();
#endif
    
    /* this mapping is only required as ISPC does not propagate task groups */

    thread2event = (ThreadEvent*) alignedMalloc(numThreads*sizeof(ThreadEvent));

    memset(thread2event,0,numThreads*sizeof(ThreadEvent));

    /* generate all threads */
    for (size_t t=0; t<numThreads; t++) {
      threads.push_back(createThread((thread_func)threadFunction,new Thread(t,numThreads,this),4*1024*1024,t));
    }

    //setAffinity(0);
    TaskLogger::init(numThreads);
  }

  void TaskScheduler::threadFunction(void* ptr) try 
  {
    Thread thread = *(Thread*) ptr;
    delete (Thread*) ptr;

    
    thread.scheduler->run(thread.threadIndex,thread.threadCount);
  }
  catch (const Terminate&) {
  }
  catch (const std::exception& e) {
    std::cout << "Error: " << e.what() << std::endl;
    exit(1);
  }

  void TaskScheduler::destroyThreads ()
  {
    terminate();
    for (size_t i=0; i<threads.size(); i++) join(threads[i]);
    threads.clear();
    alignedFree(thread2event); thread2event = NULL;
    terminateThreads = false;
  }


  // ================================================================================
  // ================================================================================
  // ================================================================================

  __aligned(64) void* volatile LockStepTaskScheduler::data = NULL;
  __aligned(64) void (* LockStepTaskScheduler::taskPtr)(void* data, const size_t threadID, const size_t numThreads) = NULL;

#if defined(__MIC__)
  __aligned(64) QuadTreeBarrier LockStepTaskScheduler::taskBarrier;
  //__aligned(64) Barrier LockStepTaskScheduler::taskBarrier;
#else
  __aligned(64) Barrier LockStepTaskScheduler::taskBarrier;
#endif

  __aligned(64) AlignedAtomicCounter32 LockStepTaskScheduler::taskCounter;

  void LockStepTaskScheduler::init(const size_t numThreads) {
    taskBarrier.init(numThreads);
  }

  void LockStepTaskScheduler::syncThreads(const size_t threadID, const size_t numThreads) {
    taskBarrier.wait(threadID,numThreads);
  }

  void LockStepTaskScheduler::syncThreadsWithReduction(const size_t threadID, 
						       const size_t numThreads,
						       void (* reductionFct)(const size_t currentThreadID,
									     const size_t childThreadID,
									     void *ptr),
						       void *ptr)
  {
    taskBarrier.syncWithReduction(threadID,numThreads,reductionFct,ptr);
  }


  void LockStepTaskScheduler::dispatchTaskMainLoop(const size_t threadID, const size_t numThreads)
  {
    while (true) {
      bool dispatch = dispatchTask(threadID,numThreads);
      if (dispatch == true) break;
    }  
  }

  bool LockStepTaskScheduler::dispatchTask(const size_t threadID, const size_t numThreads)
  {
    if (threadID == 0)
      taskCounter.reset(0);

    syncThreads(threadID, numThreads);

    if (taskPtr == NULL) 
      return true;

    (*taskPtr)((void*)data,threadID,numThreads);
    syncThreads(threadID, numThreads);
    
    return false;
  }

  void LockStepTaskScheduler::releaseThreads(const size_t numThreads)
  {
    taskPtr = NULL;
    data = NULL;
    dispatchTask(0,numThreads);  
  }

  // ================================================================================
  // ================================================================================
  // ================================================================================

  LockStepTaskScheduler4ThreadsLocalCore::LockStepTaskScheduler4ThreadsLocalCore()
  {
    taskPtr = NULL;
    data = NULL;
    for (size_t j=0;j<2;j++)
      for (size_t i=0;i<4;i++)
	threadState[j][i] = 0;
    mode = 0;
  }

  void LockStepTaskScheduler4ThreadsLocalCore::syncThreads(const size_t localThreadID) {
    const unsigned int m = mode;

    if (localThreadID == 0)
      {		
	__memory_barrier();
	threadState[m][localThreadID] = 1;
	__memory_barrier();

	while( (*(volatile unsigned int*)&threadState[m][0]) !=  0x01010101 )
	  __pause(WAIT_CYCLES);

	mode = 1 - mode;

	__memory_barrier();
	*(volatile unsigned int*)&threadState[m][0] = 0; 
      }
    else
      {
	__memory_barrier();
	threadState[m][localThreadID] = 1;
	__memory_barrier();
	
	while (threadState[m][localThreadID] == 1)
	  __pause(WAIT_CYCLES);
      }
 
  }


  void LockStepTaskScheduler4ThreadsLocalCore::dispatchTaskMainLoop(const size_t localThreadID, const size_t globalThreadID)
  {
    while (true) {
      bool dispatch = dispatchTask(localThreadID,globalThreadID);
      if (dispatch == true) break;
    }  
  }

  bool LockStepTaskScheduler4ThreadsLocalCore::dispatchTask(const size_t localThreadID, const size_t globalThreadID)
  {
    syncThreads(localThreadID);

    if (taskPtr == NULL) 
      return true;

    (*taskPtr)((void*)data,localThreadID,globalThreadID);
    syncThreads(localThreadID);
    
    return false;
  }

  void LockStepTaskScheduler4ThreadsLocalCore::releaseThreads(const size_t localThreadID, const size_t globalThreadID)
  {
    taskPtr = NULL;
    data = NULL;
    dispatchTask(localThreadID,globalThreadID);  
  }

}

