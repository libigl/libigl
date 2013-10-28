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

#ifndef __EMBREE_BARRIER_H__
#define __EMBREE_BARRIER_H__

#include "condition.h"

namespace embree
{
  /*! system barrier using operating system */
  class BarrierSys
  {
  public:

    void init(size_t count) {
      this->count = 0;
      this->full_size = count;
    }

    int wait()
    {
      count_mutex.lock();
      count++;

      if (count == full_size) {
        count = 0;
        cond.broadcast();
        count_mutex.unlock();
        return 1;
      }

      cond.wait(count_mutex);
      count_mutex.unlock();
      return 0;
    }

  protected:
    size_t count, full_size;
    MutexSys count_mutex;
    ConditionSys cond;
  };

  /* default barrier type */
  class Barrier : public BarrierSys {};
}

#endif
