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

#include "common/scene.h"
#include "primrefalloc.h"
#include "primrefblock.h"
#include "geometry/bezier1v.h"
#include "heuristic_fallback.h"

namespace embree
{
  namespace isa
  {
    /*! Generates a list of bezier curves from the scene. */
    class BezierRefGen
    {
      static const size_t maxTasks = 32;
      typedef atomic_set<PrimRefBlockT<BezierPrim> > BezierRefList;
      
    public:
      
      /*! standard constructor that schedules the task */
      BezierRefGen (size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefBlockAlloc<BezierPrim>* alloc, const Scene* scene, const size_t numTimeSteps = 1);
      
    public:
      
      /*! parallel task to iterate over the primitives */
      TASK_SET_FUNCTION(BezierRefGen,task_gen_parallel);
      
      /* input data */
    private:
      const Scene* scene;                  //!< input geometry
      const size_t numTimeSteps;
      PrimRefBlockAlloc<BezierPrim>* alloc;   //!< allocator for build primitive blocks
      
      /* intermediate data */
    private:
      TaskScheduler::Task task;
      PrimInfo pinfos[maxTasks];
      
      /* output data */
    public:
      BezierRefList prims;             //!< list of build primitives
      PrimInfo pinfo;                  //!< bounding information of primitives
    };
  }
}
