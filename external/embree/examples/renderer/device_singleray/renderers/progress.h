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

#ifndef __EMBREE_PROGRESS_H__
#define __EMBREE_PROGRESS_H__

#include "../default.h"
#include "sys/sync/atomic.h"
#include "sys/sync/mutex.h"

namespace embree
{
  /*! prints progress during rendering */
  class Progress
  {
  public:
    Progress (size_t num = 0);
    void start();
    void next();
    void end();
  private:
    void drawEmptyBar();
  private:
    MutexSys mutex;          //!< Mutex to protect progress output
    size_t curElement;       //!< Number of elements processed
    size_t numElements;      //!< Total number of elements to process
    size_t numDrawn;         //!< Number of progress characters drawn
    size_t terminalWidth;    //!< Width of terminal window in characters
  };
}

#endif
