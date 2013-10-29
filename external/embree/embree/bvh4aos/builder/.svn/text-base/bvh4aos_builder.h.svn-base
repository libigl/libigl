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

#ifndef __EMBREE_BVH4AOS_BUILDER_H__
#define __EMBREE_BVH4AOS_BUILDER_H__

#include "bvh4aos/bvh4aos.h"
#include "sys/sync/barrier.h"
#include "geometry/triangle.h"
#include "geometry/triangle_mesh.h"

namespace embree
{
  /* BVH4AOS builder. The builder is multi-threaded and implements 3
   * different build strategies: 1) Small tasks are finished in a
   * single thread (BuildTask) 2) Medium sized tasks are split into
   * two tasks using a single thread (SplitTask) and 3) Large tasks are
   * split using multiple threads on one processor. */

  class BVH4AOSBuilder : public RefCount
  {
  public:

    /*! Type of BVH built */
    typedef BVH4AOS Type;

  public:

    /*! creates the acceleration structure */
    static void create (TaskScheduler::Event* event, Accel* accel, const size_t minLeafSize = 1, const size_t maxLeafSize = inf) {
      new BVH4AOSBuilder(event,(BVH4AOS*)accel,minLeafSize,maxLeafSize,true);
    }

    static void create_fast(TaskScheduler::Event* event, Accel* accel, const size_t minLeafSize = 1, const size_t maxLeafSize = inf) {
      new BVH4AOSBuilder(event,(BVH4AOS*)accel,minLeafSize,maxLeafSize,false);
    }

    /*! Constructor. */
    BVH4AOSBuilder(TaskScheduler::Event* event, BVH4AOS* bvh, const size_t minLeafSize = 1, const size_t maxLeafSize = inf,
		   const bool highQuality = true);

  private:

    static unsigned int buildMode;

    /*! run function for BVH4AOSBuilder worker threads */
    TASK_RUN_FUNCTION(BVH4AOSBuilder,task_build_parallel);

    /*! finishes the build */
    TASK_COMPLETE_FUNCTION(BVH4AOSBuilder,finish);

  public:
    TaskScheduler::Task task;
    TriangleMesh* mesh;
    BVH4AOS* bvh; //!< Output BVH
  };
}

#endif
