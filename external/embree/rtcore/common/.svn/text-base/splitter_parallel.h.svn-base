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

#ifndef __EMBREE_SPLITTER_PARALLEL_H__
#define __EMBREE_SPLITTER_PARALLEL_H__

#include "sys/taskscheduler.h"

#include "accel.h"
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
    MultiThreadedSplitter (const TaskScheduler::ThreadInfo& thread, 
                           PrimRefAlloc* alloc, const BuildTriangle* triangles, const Vec3fa* vertices,
                           PrimRefBlockList& prims, const PrimInfo& pinfo, const Split& split, 
                           TaskScheduler::completeFunction cfun, void* cptr);
    
  public:

    /*! parallel execution of the split */
    void task_split_parallel(const TaskScheduler::ThreadInfo& thread, size_t elt); 
    static void _task_split_parallel(const TaskScheduler::ThreadInfo& thread, MultiThreadedSplitter* This, size_t elt) { This->task_split_parallel(thread,elt); }

    /*! parallel execution of spatial splits */
    void task_split_parallel_spatial(const TaskScheduler::ThreadInfo& thread, size_t elt); 
    static void _task_split_parallel_spatial(const TaskScheduler::ThreadInfo& thread, MultiThreadedSplitter* This, size_t elt) { This->task_split_parallel_spatial(thread,elt); }

    /*! reduces the parallel collected data to end up with a single best split for the left and right geometry */
    void task_split_parallel_reduce(const TaskScheduler::ThreadInfo& thread); 
    static void _task_split_parallel_reduce(const TaskScheduler::ThreadInfo& thread, MultiThreadedSplitter* This) { This->task_split_parallel_reduce(thread); }

    /* input data */
  private:
    PrimRefAlloc* alloc;             //!< allocator for build primitive blocks
    PrimRefBlockList prims;          //!< input list of build primitives
    PrimInfo pinfo;                  //!< bounding information of build primitives
    Split split;                     //!< split to perform

    /* input data for spatial splits */
  private:
    const BuildTriangle* triangles; //!< input triangle array
    const Vec3fa* vertices;         //!< input vertices
    
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
    TaskScheduler::completeFunction cfun;  //!< continuation function
    void* cptr;                            //!< call data for continuation function
  };
}

#endif

