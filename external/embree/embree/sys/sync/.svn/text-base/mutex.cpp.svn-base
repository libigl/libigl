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
}
#endif

#if defined(__UNIX__) || defined(PTHREADS_WIN32)
#include <pthread.h>
namespace embree
{
  /*! system mutex using pthreads */
  MutexSys::MutexSys( void ) 
  { 
    mutex = new pthread_mutex_t; 
    if (pthread_mutex_init((pthread_mutex_t*)mutex, NULL) != 0)
      throw std::runtime_error("pthread_mutex_init failed");
  }
  
  MutexSys::~MutexSys( void ) 
  { 
    if (pthread_mutex_destroy((pthread_mutex_t*)mutex) != 0)
      throw std::runtime_error("pthread_mutex_destroy failed");
    
    delete (pthread_mutex_t*)mutex; 
  }
  
  void MutexSys::lock( void ) 
  { 
    if (pthread_mutex_lock((pthread_mutex_t*)mutex) != 0) 
      throw std::runtime_error("pthread_mutex_lock failed");
  }
  
  void MutexSys::unlock( void ) 
  { 
    if (pthread_mutex_unlock((pthread_mutex_t*)mutex) != 0)
      throw std::runtime_error("pthread_mutex_unlock failed");
  }
}
#endif

namespace embree
{
  void MutexActive::lock () 
  {
    while (1) {
      while (flag == 1) __pause(1023); // read without atomic op first
      if (cmpxchg(flag, 1, 0) == 0) break;
    }
    __memory_barrier();  // compiler must not schedule loads and stores around this point
  }

  void MutexActive::unlock () { 
    __memory_barrier();  // compiler must not schedule loads and stores around this point
    flag = 0; 
  }
}
