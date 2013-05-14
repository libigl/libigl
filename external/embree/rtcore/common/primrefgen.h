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

#ifndef __EMBREE_PRIMREFGENBIN_H__
#define __EMBREE_PRIMREFGENBIN_H__

#include "accel.h"
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
    PrimRefGen (const TaskScheduler::ThreadInfo& thread,
                const BuildTriangle* triangles, size_t numTriangles, 
                const Vec3fa* vertices, size_t numVertices, 
                const BBox3f& bounds, PrimRefAlloc* alloc);

    /*! destruction */
    ~PrimRefGen();
    
  public:

    /*! creates a build primitive from a triangle */
    PrimRef MakePrimRef(size_t i, BBox3f& geomBound, BBox3f& centBound);

    /*! parallel task to iterate over the triangles */
    void task_gen_parallel(const TaskScheduler::ThreadInfo& thread, size_t elt); 
    static void _task_gen_parallel(const TaskScheduler::ThreadInfo& thread, PrimRefGen* This, size_t elt) { This->task_gen_parallel(thread,elt); }

    /*! reduces the bounding information gathered in parallel */
    void task_gen_parallel_reduce(const TaskScheduler::ThreadInfo& thread); 
    static void _task_gen_parallel_reduce(const TaskScheduler::ThreadInfo& thread, PrimRefGen* This) { This->task_gen_parallel_reduce(thread); }
    
    /* input data */
  private:
    const BuildTriangle* triangles;  //!< array of input triangles
    size_t numTriangles;             //!< number of triangles
    const Vec3fa* vertices;          //!< array of input vertices
    size_t numVertices;              //!< number of vertices
    PrimRefAlloc* alloc;             //!< allocator for build primitive blocks
    
    /* intermediate data */
  private:
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
