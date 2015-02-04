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
#include "sys/intrinsics.h"

namespace embree
{
#if defined(RTCORE_SPINLOCKS)
#define Barrier LinearBarrierActive
  //#define Barrier BarrierActive
#else
#define Barrier BarrierSys
#endif

  /*! system barrier using operating system */
  class BarrierSys
  {
  public:
    BarrierSys ();
    ~BarrierSys ();

  public:
    void init(size_t count);
    void wait();

  public:

    void wait (const unsigned int threadIndex, const unsigned int threadCount) {  
      wait();   
    }

    void syncWithReduction(const size_t threadIndex, 
                           const size_t threadCount,
                           void (* reductionFct)(const size_t currentThreadID,
                                                 const size_t childThreadID,
                                                 void *ptr),
                           void *ptr)
    {
      wait();
    }

  protected:
    void* opaque;
  };

  // =================================================
  // === fast memory barrier ===
  // =================================================

  struct __aligned(64) BarrierActive 
  {
  public:
    BarrierActive () 
      : cntr(0) {}
    
    void reset() {
      cntr = 0;
    }

    void wait (size_t numThreads) {
      atomic_add((atomic_t*)&cntr,1);
      while (cntr != numThreads) __pause_cpu();
    }

  private:
    volatile atomic_t cntr;
    char align[64-sizeof(atomic_t)];
  };

  // =================================================
  // === fast memory barrier with linear reduction ===
  // =================================================

#define MAX_MIC_THREADS 256 // FIXME
#define MAX_MIC_CORES (MAX_MIC_THREADS/4) // FIXME

#define MIN_MIC_BARRIER_WAIT_CYCLES 8
#define MAX_MIC_BARRIER_WAIT_CYCLES 256


  class __aligned(64) LinearBarrierActive
  {
    public:
      volatile unsigned char count0[MAX_MIC_THREADS]; 
      volatile unsigned char count1[MAX_MIC_THREADS]; 

      volatile unsigned int mode;
      volatile unsigned int flag0;
      volatile unsigned int flag1;
      volatile unsigned int fill[16-3];
      volatile unsigned int numThreads;
      
      LinearBarrierActive (size_t numThreads = 0);

      void init(size_t numThreads);

      __forceinline void pause(unsigned int &cycles) {
	__pause_cpu_expfalloff(cycles,MAX_MIC_BARRIER_WAIT_CYCLES);
      }

      void wait (const size_t threadIndex, const size_t threadCount); // FIXME: remove second parameter
      void waitForThreads(const size_t threadIndex, const size_t threadCount);

      void syncWithReduction(const size_t threadIndex, 
                             const size_t threadCount,
                             void (* reductionFct)(const size_t currentThreadID,
                                                   const size_t childThreadID,
                                                   void *ptr),
                             void *ptr);
  };



  class __aligned(64) QuadTreeBarrier
  {
  public:

    class __aligned(64) CoreSyncData {
    public:
      volatile unsigned char threadState[2][4];
      volatile unsigned int mode;
      volatile unsigned int data[16-3];

      void init();

      void pause(unsigned int &cycles);

      __forceinline void prefetchEx() { 
	prefetchL1EX((char*)threadState); 
      }

      __forceinline void prefetch() { 
	prefetchL1((char*)threadState); 
      }
    
      void switchModeAndSendRunSignal(const unsigned int m);

      void setThreadStateToDone(const unsigned int m, const unsigned int threadID);

      bool allThreadsDone(const unsigned int m, const unsigned int orMask = 0);

      bool threadDone(const unsigned int m, const unsigned int threadID);
 

      void waitForAllThreadsOnCore(const unsigned int m);

      void waitForAllOtherThreadsOnCore(const unsigned int m, const unsigned int threadID);

      void waitForThreadReceivesRunSignal(const unsigned int m, const unsigned int threadID);
    };

    QuadTreeBarrier();

    void init(size_t cntr);

    void wait(const size_t threadID, const size_t MAX_THREADS_SYNC);

    void syncWithReduction(const size_t threadID, 
                           const size_t MAX_THREADS_SYNC,
                           void (* reductionFct)(const size_t currentThreadID,
                                                 const size_t childThreadID,
                                                 void *ptr),
                           void *ptr);
    

  public:  
    __aligned(64) CoreSyncData data[MAX_MIC_CORES]; // == one cacheline per core ==
  };
}
