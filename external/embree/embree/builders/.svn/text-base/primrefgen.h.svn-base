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

#ifndef __EMBREE_PRIMREFGENBIN_H__
#define __EMBREE_PRIMREFGENBIN_H__

#include "../geometry/geometry.h"
#include "primrefalloc.h"
#include "primrefblock.h"

namespace embree
{
  /*! Generates a list of build primitives from a list of triangles. */
  template<typename Heuristic, typename PrimRefBlockList>      
    class PrimRefGen
  {
    static const size_t numTasks = 40;
    typedef typename Heuristic::Split Split;
    typedef typename Heuristic::PrimInfo PrimInfo;
    
  public:
    __forceinline PrimRefGen () {}

    /*! standard constructor that schedules the task */
    PrimRefGen (TaskScheduler::Event* group, const RTCGeometry* geom, PrimRefAlloc* alloc);

    /*! destruction */
    ~PrimRefGen();
    
  public:
    
    /*! parallel task to iterate over the triangles */
    TASK_RUN_FUNCTION(PrimRefGen,task_gen_parallel);

    /*! reduces the bounding information gathered in parallel */
    TASK_COMPLETE_FUNCTION(PrimRefGen,task_gen_parallel_reduce);
    
    /* input data */
  private:
    const RTCGeometry* geom;       //!< input geometry
    size_t numTriangles;           //!< number of primitives
    PrimRefAlloc* alloc;           //!< allocator for build primitive blocks
    
    /* intermediate data */
  private:
    TaskScheduler::Task task;
    BBox3f geomBounds[numTasks];     //!< Geometry bounds per thread
    BBox3f centBounds[numTasks];     //!< Centroid bounds per thread
    Heuristic heuristics[numTasks];  //!< Heuristics per thread
    
    /* output data */
  public:
    PrimRefBlockList prims;          //!< list of build primitives
    PrimInfo pinfo;                  //!< bounding information of primitives
    Split split;                     //!< best possible split
  };
}

#endif
