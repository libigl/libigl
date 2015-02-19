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

/*! \file Implements a task logger. One can log start and end cycle of a task and store the
   resulting scheduling diagram into a FIG file. */

#include "sys/platform.h"
#include "sys/intrinsics.h"
#include "sys/sysinfo.h"
#include <vector>

namespace embree
{
  class TaskLogger
  {
    /* maximal number of tasks and counters */
    static const size_t PERF_MAX_TASKS = 65536;

  public:
    static bool active;
    static int64 startCycle;
    static std::vector<TaskLogger*> threads;

  public:

    /* initialize the task logger */
    static bool init (size_t numThreads);

    /* start logging tasks */
    static void start();

    /* marks begin of task */
    __forceinline static size_t beginTask(size_t threadIndex, const char* name, size_t elt) 
    {
#if defined(RTCORE_TASKLOGGER)
      assert(threadIndex < threads.size());
      return threads[threadIndex]->beginTask(name,elt);
#else
      return 0;
#endif
    }

    /* marks end of task */
    __forceinline static void endTask(size_t threadIndex, size_t id) 
    {
#if defined(RTCORE_TASKLOGGER)
      assert(threadIndex < threads.size());
      threads[threadIndex]->endTask(id);
#endif
    }

    /* stops logging tasks */
    static void stop();

    /* store scheduling diagram to FIG file */
    static void store(const char* fname);

  public:
    
    TaskLogger (int threadID) : threadID(threadID) {
      reset();
    }
    
    __forceinline void reset() 
    {
      curTask = 0;
      for (size_t t=0; t<PERF_MAX_TASKS; t++) {
        counters[t].name = "uninitialized";
        counters[t].elt = 0;
        counters[t].start = 0;
        counters[t].stop = 0;
      }
    }
    
    __forceinline size_t beginTask(const char* name, size_t elt) 
    {
      if (!active || curTask >= PERF_MAX_TASKS) return PERF_MAX_TASKS-1;
      counters[curTask].name =  name;
      counters[curTask].elt = elt;
      counters[curTask].start = rdtsc() - startCycle;
      counters[curTask].stop = counters[curTask].start+1;
      return curTask++;
    }
    
    __forceinline void endTask(size_t id) 
    {
      if (!active || id >= PERF_MAX_TASKS) return;
      counters[id].stop = rdtsc() - startCycle;
    }
    
  private:
    size_t threadID;
    size_t curTask;
    struct {
      const char* name;
	  size_t elt;
      int64 start;
      int64 stop;
    } counters[PERF_MAX_TASKS];
  };
}
