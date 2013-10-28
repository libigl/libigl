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

#ifndef __EMBREE_SPLITTER_PARALLEL_H__
#define __EMBREE_SPLITTER_PARALLEL_H__

#include "sys/taskscheduler.h"
#include "../geometry/geometry.h"
#include "primrefalloc.h"
#include "primrefblock.h"

namespace embree
{
  template<typename Heuristic, typename PrimRefBlockList>      
    class MultiThreadedSplitter
  {
    static const size_t numTasks = 40;
    typedef typename Heuristic::Split Split;
    typedef typename Heuristic::PrimInfo PrimInfo;
    
  public:
    MultiThreadedSplitter () {}
    
    /*! the constructor schedules the split task */
    MultiThreadedSplitter (size_t threadIndex, size_t threadCount, TaskScheduler::Event* group, 
                           PrimRefAlloc* alloc, const RTCGeometry* geom,
                           PrimRefBlockList& prims, const PrimInfo& pinfo, const Split& split, 
                           TaskScheduler::completeFunction cfun, void* cptr);
    
  public:

    /*! parallel execution of the split */
    TASK_RUN_FUNCTION(MultiThreadedSplitter,task_split_parallel);
    
    /*! parallel execution of spatial splits */
    TASK_RUN_FUNCTION(MultiThreadedSplitter,task_split_parallel_spatial);

    /*! reduces the parallel collected data to end up with a single best split for the left and right geometry */
    TASK_COMPLETE_FUNCTION(MultiThreadedSplitter,task_split_parallel_reduce);

    /* input data */
  private:
    PrimRefAlloc* alloc;             //!< allocator for build primitive blocks
    PrimRefBlockList prims;          //!< input list of build primitives
    PrimInfo pinfo;                  //!< bounding information of build primitives
    Split split;                     //!< split to perform

    /* input data for spatial splits */
  private:
    const RTCGeometry* geom;      //!< input geometry
    
    /* intermediate data */
  private:
    Heuristic lheuristics[numTasks]; //!< parallel heuristic gathering
    Heuristic rheuristics[numTasks]; //!< parallel heuristic gathering
    
    /* output data */
  public:
    PrimRefBlockList lprims;          //!< left primitive list
    PrimRefBlockList rprims;          //!< right primitive list
    PrimInfo linfo;                   //!< left bounding information of primitives
    PrimInfo rinfo;                   //!< right bounding information of primitives
    Split lsplit;                     //!< best split for left primitives
    Split rsplit;                     //!< best split for right primitives
    
    /* continuation */
  private:
    TaskScheduler::Task task;
    TaskScheduler::completeFunction cfun;  //!< continuation function
    void* cptr;                            //!< call data for continuation function
  };
}

#endif

