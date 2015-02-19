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

#include "taskscheduler_sys.h"
#include "tasklogger.h"
 
#if 0
 
namespace embree
{
  TaskSchedulerSys::TaskSchedulerSys()
    : begin(0), end(0), tasks(16*1024) {}

  void TaskSchedulerSys::add(ssize_t threadIndex, QUEUE queue, Task* task)
  {
    if (task->event) 
      task->event->inc();

    mutex.lock();

    /*! insert task to correct end of list */
    switch (queue) {
    case GLOBAL_FRONT: { size_t i = (--begin)&(tasks.size()-1); tasks[i] = task; break; }
    case GLOBAL_BACK : { size_t i = (end++  )&(tasks.size()-1); tasks[i] = task; break; }
    default          : THROW_RUNTIME_ERROR("invalid task queue");
    }
    
    condition.broadcast();
    mutex.unlock();
  }

  void TaskSchedulerSys::run(size_t threadIndex, size_t threadCount)
  {
    while (true)
    {      
      /* wait for available task */
      mutex.lock();
      while ((end-begin) == 0 && !terminateThreads) {
        condition.wait(mutex); continue;
      }
      
      /* terminate this thread */
      if (terminateThreads) {
        mutex.unlock(); return;
      }

      /* take next task from stack */
      size_t i = (end-1)&(tasks.size()-1);
      Task* task = tasks[i]; 
      if (task == NULL || task->started <= 0) {
        tasks[i] = NULL; end--; mutex.unlock(); continue;
      }
      task->locks++;
      mutex.unlock();

      TaskScheduler::Event* event = task->event;
      thread2event[threadIndex].event = event; 

      /* run the task */
      while (true) {
        ssize_t elt = --task->started;
        if (elt < 0) break;
        if (task->run) task->run(task->runData,threadIndex,threadCount,elt,task->elts,task->event);
      }

      /* complete the task */
      mutex.lock();
      if (--task->locks == 0) {
        tasks[i] = NULL;
        mutex.unlock();
        if (task->complete) task->complete(task->completeData,threadIndex,threadCount,task->event);
        if (event) event->dec();
      } 
      else 
        mutex.unlock();
    }
  }

  void TaskSchedulerSys::terminate() 
  {
    mutex.lock();
    terminateThreads = true;
    condition.broadcast(); 
    mutex.unlock();
  }
}

#else

namespace embree
{
  TaskSchedulerSys::TaskSchedulerSys()
    : begin(0), end(0), tasks(16*1024) {}

  void TaskSchedulerSys::add(ssize_t threadIndex, QUEUE queue, Task* task)
  {
    if (task->event) 
      task->event->inc();

    mutex.lock();

    /*! resize array if too small */
    if (end-begin == tasks.size())
    {
      size_t s0 = 1*tasks.size();
      size_t s1 = 2*tasks.size();
      tasks.resize(s1);
      for (size_t i=begin; i!=end; i++)
        tasks[i&(s1-1)] = tasks[i&(s0-1)];
    }

    /*! insert task to correct end of list */
    switch (queue) {
    case GLOBAL_FRONT: { size_t i = (--begin)&(tasks.size()-1); tasks[i] = task; break; }
    case GLOBAL_BACK : { size_t i = (end++  )&(tasks.size()-1); tasks[i] = task; break; }
    default          : THROW_RUNTIME_ERROR("invalid task queue");
    }
    
    condition.broadcast();
    mutex.unlock();
  }

  void TaskSchedulerSys::wait(size_t threadIndex, size_t threadCount, Event* event)
  {
    event->dec();
    if (!isEnabled(threadIndex)) {
      while (!event->triggered()); // FIXME: wait using events
    }
    else {
      while (!event->triggered()) { // FIXME: wait using events
	work(threadIndex,threadCount,false);
      }
    }
  }

  void TaskSchedulerSys::work(size_t threadIndex, size_t threadCount, bool wait)
  {
    /* wait for available task */
    mutex.lock();
    while (((end-begin) == 0 && !terminateThreads) || !isEnabled(threadIndex)) {
      if (wait) condition.wait(mutex);
      else { mutex.unlock(); return; }
    }
    
    /* terminate this thread */
    if (terminateThreads) {
      mutex.unlock();
      return;
    }
    
    /* take next task from stack */
    size_t i = (end-1)&(tasks.size()-1);
    Task* task = tasks[i]; 
    size_t elt = --task->started;
    if (elt == 0) end--;
    mutex.unlock();
    
    /* run the task */
    TaskScheduler::Event* event = task->event;
    if (task->run) {
      size_t taskID = TaskLogger::beginTask(threadIndex,task->name,elt);
      task->run(task->runData,threadIndex,numEnabledThreads,elt,task->elts,task->event);
      TaskLogger::endTask(threadIndex,taskID);
    }
    
    /* complete the task */
    if (--task->completed == 0) {
      if (task->complete) {
        size_t taskID = TaskLogger::beginTask(threadIndex,task->name,0);
        task->complete(task->completeData,threadIndex,numEnabledThreads,task->event);
        TaskLogger::endTask(threadIndex,taskID);
      }
      if (event) event->dec();
    }
  }

  void TaskSchedulerSys::run(size_t threadIndex, size_t threadCount)
  {
    while (!terminateThreads)
      work(threadIndex,threadCount,true);
  }

  void TaskSchedulerSys::terminate() 
  {
    mutex.lock();
    terminateThreads = true;
    condition.broadcast(); 
    mutex.unlock();
  }
}

#endif
