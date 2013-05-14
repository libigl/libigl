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

#include "mutex.h"

#if defined(__WIN32__) && !defined(PTHREADS_WIN32)

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

namespace embree
{
  /*! system mutex using windows API */
  MutexSys::MutexSys( void ) { mutex = new CRITICAL_SECTION; InitializeCriticalSection((CRITICAL_SECTION*)mutex); }
  MutexSys::~MutexSys( void ) { DeleteCriticalSection((CRITICAL_SECTION*)mutex); delete (CRITICAL_SECTION*)mutex; }
  void MutexSys::lock( void ) { EnterCriticalSection((CRITICAL_SECTION*)mutex); }
  void MutexSys::unlock( void ) { LeaveCriticalSection((CRITICAL_SECTION*)mutex); }

  void MutexActive::lock  ( void ) { 
    while ( cmpxchg($lock, LOCK_IS_TAKEN, LOCK_IS_FREE) != LOCK_IS_FREE) _mm_pause();       
    __memory_barrier();  // compiler must not schedule loads and stores around this point
  }
}
#endif

#if defined(__UNIX__) || defined(PTHREADS_WIN32)
#include <pthread.h>
#include <xmmintrin.h>
namespace embree
{
  /*! system mutex using pthreads */
  MutexSys::MutexSys( void ) { mutex = new pthread_mutex_t; pthread_mutex_init((pthread_mutex_t*)mutex, NULL); }
  MutexSys::~MutexSys( void ) { pthread_mutex_destroy((pthread_mutex_t*)mutex); delete (pthread_mutex_t*)mutex; }
  void MutexSys::lock( void ) { pthread_mutex_lock((pthread_mutex_t*)mutex); }
  void MutexSys::unlock( void ) { pthread_mutex_unlock((pthread_mutex_t*)mutex); }

  void MutexActive::lock  ( void ) { 
    while ( cmpxchg($lock, LOCK_IS_TAKEN, LOCK_IS_FREE) != LOCK_IS_FREE) _mm_pause(); 
    __memory_barrier();  // compiler must not schedule loads and stores around this point
  }
}
#endif
