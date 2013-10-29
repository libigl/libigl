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
#include "taskscheduler_mic.h"
#include "sysinfo.h"

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

#if defined(__MIC__)
#if 1
    /* enable fast spinning tasking system */
    instance = new TaskSchedulerMIC;
#else
    /* enable slower pthreads tasking system */
    printf("WARNING: taskscheduler.cpp: Using pthreads tasking system active on MIC! Expect reduced rendering performance.\n");
    instance = new TaskSchedulerSys; 
#endif
#else
    /* enable fast pthreads tasking system */
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

    return instance->thread2event[threadIndex];
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
    thread2event = new Event*[numThreads];
    memset(thread2event,0,numThreads*sizeof(Event*));

    /* generate all threads */
    for (size_t t=0; t<numThreads; t++) {
      threads.push_back(createThread((thread_func)threadFunction,new Thread(t,numThreads,this),4*1024*1024,t));
    }
  }

  void TaskScheduler::threadFunction(void* ptr) try 
  {
    Thread* thread = (Thread*) ptr;
    thread->scheduler->run(thread->threadIndex,thread->threadCount);
    delete thread;
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
    delete[] thread2event; thread2event = NULL;
    terminateThreads = false;
  }
}

