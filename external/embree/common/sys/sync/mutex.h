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

#include "../platform.h"
#include "atomic.h"

namespace embree
{
  /*! system mutex */
  class MutexSys {
    friend struct ConditionImplementation;
  public:
    MutexSys( void );
    ~MutexSys( void );

    void lock( void );
    void unlock( void );

  protected:
    void* mutex;
  };

  /*! spinning mutex */
  class AtomicMutex
  {
  public:
 
    static const size_t MAX_MIC_WAIT_CYCLES = 256;
    AtomicMutex ()
      : flag(0) {}

    __forceinline bool islocked() {
      return flag == 1;
    }

    __forceinline void lock()
    {
      unsigned int wait = 128;
      while(1) {
        __memory_barrier();
	while (flag == 1) { // read without atomic op first
#if !defined(__MIC__)
	  _mm_pause(); 
	  _mm_pause();
#else
	  _mm_delay_32(wait); 
	  wait += wait;  
	  if (wait > MAX_MIC_WAIT_CYCLES) wait = MAX_MIC_WAIT_CYCLES;  
#endif
	}
        __memory_barrier();
	if (atomic_cmpxchg(&flag,0,1) == 0) break;
#if defined(__MIC__)
	_mm_delay_32(wait); 
#endif
      }
    }

    __forceinline void unlock() 
    {
      __memory_barrier();
      flag = 0;
      __memory_barrier();
    }

    __forceinline void reset(int i = 0) 
    {
      __memory_barrier();
      flag = i;
      __memory_barrier();
    }

  public:
    volatile int flag;
  };

  class __aligned(64) AlignedAtomicMutex : public AtomicMutex
  {
  public:
    volatile unsigned int index;
    atomic32_t m_counter;
    volatile char align[64-3*sizeof(int)]; 
  
    AlignedAtomicMutex() {
      m_counter = 0;
      index = 0;
    }

    __forceinline void resetCounter(unsigned int i = 0) {
        __memory_barrier();
      *(volatile unsigned int*)&m_counter = i;
        __memory_barrier();
    }

    __forceinline unsigned int inc() {
      return atomic_add(&m_counter,1);
    }

    __forceinline unsigned int dec() {
      return atomic_add(&m_counter,-1);
    }

    __forceinline unsigned int val() {
      __memory_barrier();
      return m_counter;
    };
  };

  /*! safe mutex lock and unlock helper */
  template<typename Mutex> class Lock {
  public:
    Lock (Mutex& mutex) : mutex(mutex) { mutex.lock(); }
    ~Lock() { mutex.unlock(); }
  protected:
    Mutex& mutex;
  };
}
