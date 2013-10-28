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

#include "condition.h"

#if defined(__WIN32__) && !defined(PTHREADS_WIN32)

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

namespace embree
{
  /* 
     This is an implementation of POSIX "compatible" condition variables 
     for Win32, as described by Douglas C. Schmidt and Irfan Pyarali: 
     www.cs.wustl.edu/~schmidt/win32-cv-1.html. The code is a modified 
     version from the glfw source code (www.glfw.org) published under 
     the zlib/libpng license:

     Copyright © 2002-2006 Marcus Geelnard
     Copyright © 2006-2012 Camilla Berglund
     
     This software is provided 'as-is', without any express or implied
     warranty. In no event will the authors be held liable for any
     damages arising from the use of this software.
     
     Permission is granted to anyone to use this software for any
     purpose, including commercial applications, and to alter it and
     redistribute it freely, subject to the following restrictions:
     
     1. The origin of this software must not be misrepresented; you
     must not claim that you wrote the original software. If you use
     this software in a product, an acknowledgment in the product
     documentation would be appreciated but is not required.
     
     2. Altered source versions must be plainly marked as such, and
     must not be misrepresented as being the original software.
     
     3. This notice may not be removed or altered from any source
     distribution.
  */

  enum {
      COND_SIGNAL     = 0,
      COND_BROADCAST  = 1
  };

  /*! Implement the internal condition variable implementation Mingw lacks */
  struct ConditionWindows
  {
    HANDLE events[2];                     //<! Signal and broadcast event HANDLEs
    unsigned int waiters_count;           //<! Count of the number of waiters
    CRITICAL_SECTION waiters_count_lock;  //!< Serialize access to <waiters_count>
  };

  ConditionSys::ConditionSys ()
  {
    cond = new ConditionWindows;
    ConditionWindows* cv = (ConditionWindows*) cond;
    cv->waiters_count = 0;
    cv->events[COND_SIGNAL]    = CreateEvent(NULL, FALSE, FALSE, NULL);
    cv->events[COND_BROADCAST] = CreateEvent(NULL, TRUE, FALSE, NULL);
    InitializeCriticalSection(&cv->waiters_count_lock);
  }

  ConditionSys::~ConditionSys ()
  {
    ConditionWindows* cv = (ConditionWindows*) cond;
    CloseHandle(cv->events[COND_SIGNAL]);
    CloseHandle(cv->events[COND_BROADCAST]);
    DeleteCriticalSection(&cv->waiters_count_lock);
    delete cv;
  }

  void ConditionSys::wait(MutexSys& mutex)
  {
    ConditionWindows* cv = (ConditionWindows*) cond;
    int result, last_waiter;
    DWORD timeout_ms;

    // Avoid race conditions
    EnterCriticalSection(&cv->waiters_count_lock);
    cv->waiters_count ++;
    LeaveCriticalSection(&cv->waiters_count_lock);

    // It's ok to release the mutex here since Win32 manual-reset events
    // maintain state when used with SetEvent()
    LeaveCriticalSection((CRITICAL_SECTION *) mutex.mutex);
    timeout_ms = INFINITE;

    // Wait for either event to become signaled due to glfwSignalCond or
    // glfwBroadcastCond being called
    result = WaitForMultipleObjects(2, cv->events, FALSE, timeout_ms);

    // Check if we are the last waiter
    EnterCriticalSection(&cv->waiters_count_lock);
    cv->waiters_count --;
    last_waiter = (result == WAIT_OBJECT_0 + COND_BROADCAST) &&
                  (cv->waiters_count == 0);
    LeaveCriticalSection(&cv->waiters_count_lock);

    // Some thread called glfwBroadcastCond
    if (last_waiter) {
      // We're the last waiter to be notified or to stop waiting, so
      // reset the manual event
      ResetEvent(cv->events[COND_BROADCAST]);
    }

    // Reacquire the mutex
    EnterCriticalSection((CRITICAL_SECTION *) mutex.mutex);
  }

  void ConditionSys::broadcast()
  {
    ConditionWindows* cv = (ConditionWindows*) cond;
    int have_waiters;

    // Avoid race conditions
    EnterCriticalSection(&cv->waiters_count_lock);
    have_waiters = cv->waiters_count > 0;
    LeaveCriticalSection(&cv->waiters_count_lock);

    if (have_waiters)
      SetEvent(cv->events[COND_BROADCAST]);
  }
}
#endif

#if defined(__UNIX__) || defined(PTHREADS_WIN32)
#include <pthread.h>
namespace embree
{
  ConditionSys::ConditionSys () { cond = new pthread_cond_t; pthread_cond_init((pthread_cond_t*)cond,NULL); }
  ConditionSys::~ConditionSys() { pthread_cond_destroy((pthread_cond_t*)cond); delete (pthread_cond_t*)cond; } 
  void ConditionSys::wait(MutexSys& mutex) { pthread_cond_wait((pthread_cond_t*)cond, (pthread_mutex_t*)mutex.mutex); }
  void ConditionSys::broadcast() { pthread_cond_broadcast((pthread_cond_t*)cond); }
}
#endif
