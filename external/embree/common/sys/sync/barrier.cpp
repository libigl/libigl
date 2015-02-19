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

#include "barrier.h"
#include "condition.h"

#if defined (__WIN32__)

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

namespace embree
{
  struct BarrierSysImplementation
  {
    __forceinline BarrierSysImplementation () 
      : i(0), enterCount(0), exitCount(0), barrierSize(0) 
    {
      events[0] = CreateEvent(NULL, TRUE, FALSE, NULL);
      events[1] = CreateEvent(NULL, TRUE, FALSE, NULL);
    }
    
    __forceinline ~BarrierSysImplementation ()
    {
      CloseHandle(events[0]);
      CloseHandle(events[1]);
    }
    
    __forceinline void init(size_t N) 
    {
      barrierSize = N;
      enterCount = N;
      exitCount = N;
    }

    __forceinline void wait()
    {
      /* every thread entering the barrier decrements this count */
      size_t i0 = i;
      ssize_t cnt0 = atomic_add(&enterCount,-1);

      /* all threads except the last one are wait in the barrier */
      if (cnt0 > 1) 
      {
        if (WaitForSingleObject(events[i0], INFINITE) != WAIT_OBJECT_0)
          THROW_RUNTIME_ERROR("WaitForSingleObjects failed");
      }
      
      /* the last thread starts all threads waiting at the barrier */
      else 
      {
        i = 1-i;
        enterCount = barrierSize;
        if (SetEvent(events[i0]) == 0)
          THROW_RUNTIME_ERROR("SetEvent failed");
      }

      /* every thread leaving the barrier decrements this count */
      ssize_t cnt1 = atomic_add(&exitCount,-1);

      /* the last thread that left the barrier resets the event again */
      if (cnt1 == 1) 
      {
        exitCount = barrierSize;
        if (ResetEvent(events[i0]) == 0)
          THROW_RUNTIME_ERROR("ResetEvent failed");
      }
    }

  public:
    HANDLE events[2];
    volatile size_t i;
    volatile atomic_t enterCount;
    volatile atomic_t exitCount;
	size_t barrierSize;
  };
}

#else

namespace embree
{
  struct BarrierSysImplementation
  {
    __forceinline BarrierSysImplementation () 
      : count(0), barrierSize(0) {}
    
    __forceinline void init(size_t N) 
    {
      count = 0;
      barrierSize = N;
    }

    __forceinline void wait()
    {
      mutex.lock();
      count++;
      
      if (count == barrierSize) {
        count = 0;
        cond.broadcast();
        mutex.unlock();
        return;
      }
      
      cond.wait(mutex);
      mutex.unlock();
      return;
    }

  public:
    MutexSys mutex;
    ConditionSys cond;
    volatile size_t count;
    volatile size_t barrierSize;
  };
}

#endif

namespace embree
{
  BarrierSys::BarrierSys () {
    opaque = new BarrierSysImplementation;
  }

  BarrierSys::~BarrierSys () {
    delete (BarrierSysImplementation*) opaque;
  }

  void BarrierSys::init(size_t count) {
    ((BarrierSysImplementation*) opaque)->init(count);
  }

  void BarrierSys::wait() {
    ((BarrierSysImplementation*) opaque)->wait();
  }

  LinearBarrierActive::LinearBarrierActive (size_t numThreads_i) 
  {    
    numThreads = numThreads_i;
    mode      = 0;
    flag0     = 0;
    flag1     = 0;
    for (size_t i=0; i<MAX_MIC_THREADS; i++) count0[i] = 0;
    for (size_t i=0; i<MAX_MIC_THREADS; i++) count1[i] = 0;
  }

  void LinearBarrierActive::init(size_t cntr) 
  {
    numThreads = cntr;
    mode      = 0;
    flag0     = 0;
    flag1     = 0;
    for (size_t i=0; i<cntr; i++) count0[i] = 0;
    for (size_t i=0; i<cntr; i++) count1[i] = 0;
  }

  void LinearBarrierActive::wait (const size_t threadIndex, const size_t __threadCount)
  {
    waitForThreads(threadIndex,numThreads);
  }

  void LinearBarrierActive::waitForThreads(const size_t threadIndex, const size_t threadCount)
  {
    if (mode == 0)
    {			
      if (threadIndex == 0)
      {	
        for (size_t i=0; i<threadCount; i++)
          count1[i] = 0;
        
        for (size_t i=1; i<threadCount; i++)
        {
          unsigned int wait_cycles = MIN_MIC_BARRIER_WAIT_CYCLES;
          while (likely(count0[i] == 0)) 
          {
            pause(wait_cycles);
          }
        }
        mode  = 1;
        flag1 = 0;
        __memory_barrier();
        flag0 = 1;
      }			
      else
      {					
        count0[threadIndex] = 1;
        {
          unsigned int wait_cycles = MIN_MIC_BARRIER_WAIT_CYCLES;
          while (likely(flag0 == 0))
          {
            pause(wait_cycles);
          }
        }
        
      }		
    }					
    else						
    {
      if (threadIndex == 0)
      {	
        for (size_t i=0; i<threadCount; i++)
          count0[i] = 0;
        
        for (size_t i=1; i<threadCount; i++)
        {		
          unsigned int wait_cycles = MIN_MIC_BARRIER_WAIT_CYCLES;
          while (likely(count1[i] == 0))
          {
            pause(wait_cycles);
          }
        }
        
        mode  = 0;
        flag0 = 0;
        __memory_barrier();
        flag1 = 1;
      }			
      else
      {					
        count1[threadIndex] = 1;
        {
          unsigned int wait_cycles = MIN_MIC_BARRIER_WAIT_CYCLES;
          while (likely(flag1 == 0))
          {
            pause(wait_cycles);
          }
        }
      }		
    }					
  }

  
  void LinearBarrierActive::syncWithReduction(const size_t threadIndex, 
                                              const size_t threadCount,
                                              void (* reductionFct)(const size_t currentThreadID,
                                                                    const size_t childThreadID,
                                                                    void *ptr),
                                              void *ptr)
  {
    if (mode == 0)
    {			
      if (threadIndex == 0)
      {	
        for (size_t i=0; i<threadCount; i++)
          count1[i] = 0;
        
        for (size_t i=1; i<threadCount; i++)
        {
          unsigned int wait_cycles = MIN_MIC_BARRIER_WAIT_CYCLES;
          while (likely(count0[i] == 0)) 
          {
            pause(wait_cycles);
          }
          (*reductionFct)(threadIndex,i,ptr);
        }
        mode  = 1;
        flag1 = 0;
        __memory_barrier();
        flag0 = 1;
      }			
      else
      {					
        count0[threadIndex] = 1;
        {
          unsigned int wait_cycles = MIN_MIC_BARRIER_WAIT_CYCLES;
          while (likely(flag0 == 0))
          {
            pause(wait_cycles);
          }
        }
      }		
    }					
    else						
    {
      if (threadIndex == 0)
      {	
        for (size_t i=0; i<threadCount; i++)
          count0[i] = 0;
        
        for (size_t i=1; i<threadCount; i++)
        {		
          unsigned int wait_cycles = MIN_MIC_BARRIER_WAIT_CYCLES;
          while (likely(count1[i] == 0))
          {
            pause(wait_cycles);
          }
          (*reductionFct)(threadIndex,i,ptr);
        }
        
        mode  = 0;
        flag0 = 0;
        __memory_barrier();
        flag1 = 1;
      }			
      else
      {					
        count1[threadIndex] = 1;
        {
          unsigned int wait_cycles = MIN_MIC_BARRIER_WAIT_CYCLES;
          while (likely(flag1 == 0))
          {
            pause(wait_cycles);
          }
        }
      }		
    }				             	
  }

  void QuadTreeBarrier::CoreSyncData::init() 
  {
    *(volatile unsigned int*)&threadState[0][0] = 0; 
    *(volatile unsigned int*)&threadState[1][0] = 0; 
    mode = 0;
    data[0] = 0;
    __memory_barrier();
  }
  
  void QuadTreeBarrier::CoreSyncData::pause(unsigned int &cycles) {
    __pause_cpu_expfalloff(cycles,MAX_MIC_BARRIER_WAIT_CYCLES);
  }
  
  void QuadTreeBarrier::CoreSyncData::switchModeAndSendRunSignal(const unsigned int m)
  { 
    //__memory_barrier();
    mode = 1 - mode;
    __memory_barrier();
    *(volatile unsigned int*)&threadState[m][0] = 0; 
    //__memory_barrier();
  }

  void QuadTreeBarrier::CoreSyncData::setThreadStateToDone(const unsigned int m, const unsigned int threadID)
  {
    __memory_barrier();
    threadState[m][threadID % 4] = 1;
    __memory_barrier();
  }

  bool QuadTreeBarrier::CoreSyncData::allThreadsDone(const unsigned int m, const unsigned int orMask) { 
    return (*(volatile unsigned int*)&threadState[m][0] | orMask)== 0x01010101; 
  }

  bool QuadTreeBarrier::CoreSyncData::threadDone(const unsigned int m, const unsigned int threadID) { 
    return threadState[m][threadID % 4] == 1; 
  }
  
  void QuadTreeBarrier::CoreSyncData::waitForAllThreadsOnCore(const unsigned int m)
  {
    unsigned int count = MIN_MIC_BARRIER_WAIT_CYCLES;
    while(likely(allThreadsDone(m) == false)) 
      pause(count);
  }

  void QuadTreeBarrier::CoreSyncData::waitForAllOtherThreadsOnCore(const unsigned int m, const unsigned int threadID)
  {
    unsigned int count = MIN_MIC_BARRIER_WAIT_CYCLES;
    const unsigned int orMask = (unsigned int)1 << ((threadID % 4) * 8);
    while(likely(allThreadsDone(m,orMask) == false)) 
      pause(count);
  }

  QuadTreeBarrier::QuadTreeBarrier()
  {    
    assert(sizeof(CoreSyncData) == 64);
    for (size_t i=0;i<MAX_MIC_CORES;i++) 
    {
      data[i].init();
    }
  }
  
  void QuadTreeBarrier::init(size_t cntr) {
  }
  
  void QuadTreeBarrier::CoreSyncData::waitForThreadReceivesRunSignal(const unsigned int m, const unsigned int threadID)
  {
    unsigned int count = MIN_MIC_BARRIER_WAIT_CYCLES;
    while(likely(threadDone(m,threadID) == true)) 
      pause(count);
  }
  
  void QuadTreeBarrier::wait(const size_t threadID, const size_t MAX_THREADS_SYNC)
  {
    if (unlikely(MAX_THREADS_SYNC == 1)) return;
    
    const unsigned int MAX_CORES_SYNC = MAX_THREADS_SYNC >> 2;
    const unsigned int coreID = threadID >> 2;
    const unsigned int MODE = data[coreID].mode;
    
    // Drain store buffer for NGO stores
    //atomic_add((atomic_t*)&data[coreID].data[0],0);
    
    //data[coreID].prefetchEx(); 
    
    if (threadID == 0)
    {		
      data[0].setThreadStateToDone(MODE,threadID);
      
      // == wait for core 0 ==
      data[0].waitForAllThreadsOnCore(MODE);
      
      // == wait for (possible) two children cores
      const unsigned int nextCoreID0 = 1;
      const unsigned int nextCoreID1 = 2;
      
      data[nextCoreID0].prefetch(); 
      data[nextCoreID1].prefetch(); 
      
      if (nextCoreID0 < MAX_CORES_SYNC) data[nextCoreID0].waitForAllThreadsOnCore(MODE);
      if (nextCoreID1 < MAX_CORES_SYNC) data[nextCoreID1].waitForAllThreadsOnCore(MODE);
      
      data[nextCoreID0].prefetchEx(); 
      data[nextCoreID1].prefetchEx(); 
      
      // == run signal to core 0 ==
      data[0].switchModeAndSendRunSignal(MODE);
      
      // == propagate run signal to core 1,2 ==
      if (nextCoreID0 < MAX_CORES_SYNC) data[nextCoreID0].switchModeAndSendRunSignal(MODE);
      if (nextCoreID1 < MAX_CORES_SYNC) data[nextCoreID1].switchModeAndSendRunSignal(MODE);
    }			
    else
    {				

      const unsigned int nextCoreID0 = 2*coreID + 1;
      const unsigned int nextCoreID1 = 2*coreID + 2;
      
      data[nextCoreID0].prefetch(); 
      data[nextCoreID1].prefetch(); 
      
      if (threadID % 4 == 0)
      {	
        if (nextCoreID0 < MAX_CORES_SYNC) data[nextCoreID0].waitForAllThreadsOnCore(MODE);
        if (nextCoreID1 < MAX_CORES_SYNC) data[nextCoreID1].waitForAllThreadsOnCore(MODE);
      }
      
      data[coreID].setThreadStateToDone(MODE,threadID % 4);
      data[coreID].waitForThreadReceivesRunSignal(MODE,threadID % 4);	
      
      // == propagte run signal to the two children ==
      
      if (threadID % 4 == 0)
      {
        data[nextCoreID0].prefetchEx(); 
        data[nextCoreID1].prefetchEx(); 
        
        if (nextCoreID0 < MAX_CORES_SYNC) data[nextCoreID0].switchModeAndSendRunSignal(MODE);
        if (nextCoreID1 < MAX_CORES_SYNC) data[nextCoreID1].switchModeAndSendRunSignal(MODE);
      }
    }
  }
  
  void QuadTreeBarrier::syncWithReduction(const size_t threadID, 
                                          const size_t MAX_THREADS_SYNC,
                                          void (* reductionFct)(const size_t currentThreadID,
                                                                const size_t childThreadID,
                                                                void *ptr),
                                          void *ptr)
  {
    if (unlikely(MAX_THREADS_SYNC == 1)) return;
    
    const unsigned int MAX_CORES_SYNC = MAX_THREADS_SYNC >> 2;
    const unsigned int coreID = threadID >> 2;
    const unsigned int MODE = data[coreID].mode;
    
    // Drain store buffer for NGO stores
    //atomic_add((atomic_t*)&data[coreID].data[0],0);
    
    data[coreID].prefetchEx(); 
    
    if (threadID == 0)
    {		
      
      data[0].setThreadStateToDone(MODE,threadID);
      
      // == wait for core 0 ==
      data[0].waitForAllThreadsOnCore(MODE);
      
      (*reductionFct)(threadID,threadID+1,ptr);
      (*reductionFct)(threadID,threadID+2,ptr);
      (*reductionFct)(threadID,threadID+3,ptr);
      
      // == wait for (possible) two children cores
      const unsigned int nextCoreID0 = 1;
      const unsigned int nextCoreID1 = 2;
      const unsigned int nextThreadID0 = nextCoreID0 * 4;
      const unsigned int nextThreadID1 = nextCoreID1 * 4;
      
      data[nextCoreID0].prefetch(); 
      data[nextCoreID1].prefetch(); 
      
      if (nextCoreID0 < MAX_CORES_SYNC) 
      {
        data[nextCoreID0].waitForAllThreadsOnCore(MODE);
        (*reductionFct)(threadID,nextThreadID0,ptr);
      }
      
      if (nextCoreID1 < MAX_CORES_SYNC) 
      {
        data[nextCoreID1].waitForAllThreadsOnCore(MODE);
        (*reductionFct)(threadID,nextThreadID1,ptr);
      }
      
      data[nextCoreID0].prefetchEx(); 
      data[nextCoreID1].prefetchEx(); 
      
      // == run signal to core 0 ==
      data[0].switchModeAndSendRunSignal(MODE);
      
      // == propagate run signal to core 1,2 ==
      if (nextCoreID0 < MAX_CORES_SYNC) data[nextCoreID0].switchModeAndSendRunSignal(MODE);
      if (nextCoreID1 < MAX_CORES_SYNC) data[nextCoreID1].switchModeAndSendRunSignal(MODE);
    }			
    else
    {				
      
      const unsigned int nextCoreID0 = 2*coreID + 1;
      const unsigned int nextCoreID1 = 2*coreID + 2;
      const unsigned int nextThreadID0 = nextCoreID0 * 4;
      const unsigned int nextThreadID1 = nextCoreID1 * 4;
      
      data[nextCoreID0].prefetch(); 
      data[nextCoreID1].prefetch(); 
      
      if (threadID % 4 == 0)
      {	
        data[coreID].waitForAllOtherThreadsOnCore(MODE,threadID);
        (*reductionFct)(threadID,threadID+1,ptr);
        (*reductionFct)(threadID,threadID+2,ptr);
        (*reductionFct)(threadID,threadID+3,ptr);
        
        if (nextCoreID0 < MAX_CORES_SYNC) 
        {
          data[nextCoreID0].waitForAllThreadsOnCore(MODE);
          (*reductionFct)(threadID,nextThreadID0,ptr);
        }
	
        
        if (nextCoreID1 < MAX_CORES_SYNC) 
        {
          data[nextCoreID1].waitForAllThreadsOnCore(MODE);
          (*reductionFct)(threadID,nextThreadID1,ptr);
        }
      }
      
      data[coreID].setThreadStateToDone(MODE,threadID % 4);
      
      data[coreID].waitForThreadReceivesRunSignal(MODE,threadID % 4);	
      
      // == propagte run signal to the two children ==
      
      if (threadID % 4 == 0)
      {
        data[nextCoreID0].prefetchEx(); 
        data[nextCoreID1].prefetchEx(); 
        
        if (nextCoreID0 < MAX_CORES_SYNC) data[nextCoreID0].switchModeAndSendRunSignal(MODE);
        if (nextCoreID1 < MAX_CORES_SYNC) data[nextCoreID1].switchModeAndSendRunSignal(MODE);
      }
    }        
  }
}
