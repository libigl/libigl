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

#include "taskscheduler_standard.h"
#include "sysinfo.h"
#include <xmmintrin.h>

namespace embree
{
  TaskScheduler* scheduler = new TaskSchedulerStandard;
  
  TaskSchedulerStandard::TaskSchedulerStandard (size_t numThreads_in) : activeTasks(0)
  {
    this->numThreads = numThreads_in;
    if (numThreads == 0) numThreads = getNumberOfLogicalThreads();
    terminateThreads = false; 

    /* generate all threads */
    for (size_t t=0; t<numThreads; t++)
      threads.push_back(createThread((thread_func)threadFunction,new Thread(t,this),4*1024*1024,t));
  }

  TaskSchedulerStandard::~TaskSchedulerStandard()
  {
    if (threads.size() == 0) return;
    terminateThreads = true;
    taskMutex.lock(); 
    taskCondition.broadcast(); 
    taskMutex.unlock();
    for (size_t i=0; i<threads.size(); i++) join(threads[i]);
    threads.clear();
    terminateThreads = false;
  }

  void TaskSchedulerStandard::addTask(const ThreadInfo& thread, QUEUE queue, 
                                      runFunction run, void* runData, size_t elts, completeFunction complete, void* completeData, const char* name)
  {
    /* create new task */
    Task* task = new Task(run,runData,elts,complete,completeData,name);
    activeTasks++;
    
    /* add new task to task stack */
    Lock<MutexSys> lock(taskMutex);
    switch (queue) {
    case GLOBAL_FRONT: tasks.push_front(task); break;
    case GLOBAL_BACK : tasks.push_back (task); break;
    default     : throw std::runtime_error("invalid task queue");
    }
    taskCondition.broadcast();
  }

  void TaskSchedulerStandard::start() {
    mainMutex.lock();
  }

  void TaskSchedulerStandard::stop() 
  {
    finishMutex.lock();
    while (activeTasks != 0) finishCondition.wait(finishMutex);
    finishMutex.unlock();
    mainMutex.unlock();
  }

  void TaskSchedulerStandard::run(size_t tid)
  {
    ThreadInfo tinfo(tid);
    
    while (true)
    {
      /* wait until a new task is available */
      taskMutex.lock();
      while (tasks.empty() && !terminateThreads) taskCondition.wait(taskMutex);
      if (terminateThreads) { taskMutex.unlock(); break; }
  
      /* take next task from stack */
      Task* task = tasks.front(); //top;
      size_t elt = --task->started;
      if (elt == 0) tasks.pop_front(); //top = task->next;
      taskMutex.unlock();
      
      /* run the task */
      if (task->run) task->run(tinfo,task->runData,elt);
      
      /* complete the task */
      if (--task->completed == 0) {
        if (task->complete) task->complete(tinfo,task->completeData);
        delete task;
        if (--activeTasks == 0) {
          finishMutex.lock();
          finishCondition.broadcast();
          finishMutex.unlock();
        }
      }
    }
  }

  void TaskSchedulerStandard::threadFunction(Thread* thread) try 
  {
    thread->scheduler->run(thread->tid);
    delete thread;
  }
  catch (const std::exception& e) {
    std::cout << "Error: " << e.what() << std::endl;
    exit(1);
  }
}

