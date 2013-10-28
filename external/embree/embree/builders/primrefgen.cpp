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

#include "primrefgen.h"
#include "heuristics.h"

namespace embree
{
  template<typename Heuristic, typename PrimRefBlockList>
  PrimRefGen<Heuristic,PrimRefBlockList>::PrimRefGen(TaskScheduler::Event* event, const RTCGeometry* geom, PrimRefAlloc* alloc)
    : geom(geom), numTriangles(geom->size()), alloc(alloc)
  {
    /* approximate bounds if not available */
    BBox3f approxBounds = geom->approx;
    if (approxBounds.empty() && numTriangles) {
      BBox3f geomBounds = empty, centBounds = empty;
      for (size_t i=0; i<2048; i++) {
        const size_t j = i*numTriangles/2048;
        const BBox3f b = geom->bounds(j);
        if (isEmpty(b)) continue;
        geomBounds.grow(b);
        centBounds.grow(center2(b));
      }
      new (&pinfo) PrimInfo(numTriangles,geomBounds,centBounds);
    }
    else
      new (&pinfo) PrimInfo(numTriangles,approxBounds);

    /* start parallel task */
    new (&task) TaskScheduler::Task(event,
                                    _task_gen_parallel,this,numTasks,
                                    _task_gen_parallel_reduce,this,
                                    "build::primrefgen");
    TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_FRONT,&task);
  }

  template<typename Heuristic, typename PrimRefBlockList>
  PrimRefGen<Heuristic,PrimRefBlockList>::~PrimRefGen() {
    pinfo.clear();
  }
    
  template<typename Heuristic, typename PrimRefBlockList>
  void PrimRefGen<Heuristic,PrimRefBlockList>::task_gen_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event) 
  {
    Heuristic& heuristic = heuristics[taskIndex];
    new (&heuristic) Heuristic(pinfo,geom);
    
    /* static work allocation */
    size_t start = (taskIndex+0)*numTriangles/taskCount;
    size_t end   = (taskIndex+1)*numTriangles/taskCount;
    BBox3f geomBound = empty, centBound = empty;
    typename PrimRefBlockList::item* block = prims.insert(alloc->malloc(threadIndex)); 

    for (size_t i=start; i<end; i++) 
    {
      const BBox3f b = geom->bounds(i);
      if (isEmpty(b)) continue;
      geomBound.grow(b);
      centBound.grow(center2(b));
      const PrimRef prim = PrimRef(b,i);
      if (likely(block->insert(prim))) continue; 
      heuristic.bin(block->base(),block->size());
      block = prims.insert(alloc->malloc(threadIndex));
      block->insert(prim);
    }
    heuristic.bin(block->base(),block->size());
    geomBounds[taskIndex] = geomBound;
    centBounds[taskIndex] = centBound;
  }
  
  template<typename Heuristic, typename PrimRefBlockList>
  void PrimRefGen<Heuristic,PrimRefBlockList>::task_gen_parallel_reduce(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event) 
  {
    /* reduce geometry and centroid bounds */
    BBox3f geomBound = empty;
    BBox3f centBound = empty;
    for (size_t i=0; i<numTasks; i++) {
      geomBound = merge(geomBound,geomBounds[i]);
      centBound = merge(centBound,centBounds[i]);
    }
    new (&pinfo) PrimInfo(numTriangles,geomBound,centBound);
    
    /* reduce heuristic and find best split */
    Heuristic heuristic; 
    Heuristic::reduce(heuristics,numTasks,heuristic); 
    heuristic.best(split);
  }
  
  /*! explicit template instantiations */
  INSTANTIATE_TEMPLATE_BY_HEURISTIC_AND_PRIMREFLIST(PrimRefGen);
}
