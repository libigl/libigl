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

#include "mutex.h"
#include "condition.h"

namespace embree
{
  
  class EventSys
  {
  public:
    EventSys() 
      : event(false) {}
    
    void reset() {
      event = false;
    }

    void signal() {
#if defined(__MIC__)
      event = true;
#else
      mutex.lock();
      event = true;
      condition.broadcast(); // this broadcast has to be protected!
      mutex.unlock();
#endif
    }

    void wait() {
#if defined(__MIC__)
      while (!event) __pause_cpu(1024);
#else
      mutex.lock();
      while (!event) condition.wait(mutex);
      mutex.unlock();
#endif

    }

  protected:
    volatile bool event;
#if !defined(__MIC__)
    MutexSys mutex;
    ConditionSys condition;
#endif
  };
}
