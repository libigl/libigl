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
  template<typename Heuristic>
  PrimRefGen<Heuristic>::PrimRefGen(size_t threadIndex, size_t threadCount, const BuildSource* geom, PrimRefAlloc* alloc)
    : geom(geom), numPrimitives(0), numVertices(0), alloc(alloc)
  {
    /* compute number of primitives */
    size_t numGroups = geom->groups();
    for (size_t g=0; g<numGroups; g++) {
      size_t vertices = 0;
      numPrimitives += geom->prims(g,&vertices);
      numVertices += vertices;
    }

    /* approximate bounds */
    BBox3f geomBound = empty, centBound = empty;
    size_t s = 0, t = 0, dt = max(size_t(1),numPrimitives/2048);
    for (size_t g=0; g<numGroups; g++) 
    {
      size_t numPrims = geom->prims(g);
      for (size_t i=t-s; i<numPrims; i+=dt, t+=dt) {
        BBox3f bounds = geom->bounds(g,i);
        geomBound.extend(bounds);
        centBound.extend(center2(bounds));
      }
      s += numPrims;
    }
    new (&pinfo) PrimInfo(numPrimitives,geomBound,centBound);

    /* compute start group and primitives */
    size_t g=0, i=0;
    for (size_t k=0; k<numTasks; k++) 
    {
      size_t start = (k+0)*numPrimitives/numTasks;
      size_t end   = (k+1)*numPrimitives/numTasks;
      size_t size  = end-start;
      work[k].startGroup = g;
      work[k].startPrim = i;
      work[k].numPrims = size;

      for (; g<numGroups; g++) 
      {
        size_t numPrims = geom->prims(g)-i;
        if (size < numPrims) {
          i += size;
          break;
        }
        size -= numPrims;
        i = 0;
      }
    }

    /* start parallel task */
    TaskScheduler::executeTask(threadIndex,threadCount,
                               _task_gen_parallel,this,numTasks,
                               _task_gen_parallel_reduce,this,
                               "build::primrefgen");
  }

  template<typename Heuristic>
  PrimRefGen<Heuristic>::~PrimRefGen() {
    pinfo.clear();
  }
    
  template<typename Heuristic>
  void PrimRefGen<Heuristic>::task_gen_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event) 
  {
    Heuristic& heuristic = heuristics[taskIndex];
    new (&heuristic) Heuristic(pinfo,geom);
    
    /* static work allocation */
    size_t g = work[taskIndex].startGroup;
    size_t i = work[taskIndex].startPrim;
    size_t numPrims = work[taskIndex].numPrims;
    size_t numGroupPrims = numPrims ? geom->prims(g) : 0;
    size_t numAddedPrims = 0;
    
    BBox3f geomBound = empty, centBound = empty;
    typename atomic_set<PrimRefBlock>::item* block = prims.insert(alloc->malloc(threadIndex)); 
    for (size_t p=0; p<numPrims; p++, i++)
    {
      /* goto next group */
      while (i == numGroupPrims) {
        g++; i = 0;
        numGroupPrims = geom->prims(g);
      }

      const BBox3f b = geom->bounds(g,i);
      if (b.empty()) continue;
      numAddedPrims++;
      geomBound.extend(b);
      centBound.extend(center2(b));
      const PrimRef prim = PrimRef(b,g,i);
      if (likely(block->insert(prim))) continue; 
      heuristic.bin(block->base(),block->size());
      block = prims.insert(alloc->malloc(threadIndex));
      block->insert(prim);
    }
    heuristic.bin(block->base(),block->size());
    geomBounds[taskIndex] = geomBound;
    centBounds[taskIndex] = centBound;
    work[taskIndex].numPrims = numAddedPrims;
  }
  
  template<typename Heuristic>
  void PrimRefGen<Heuristic>::task_gen_parallel_reduce(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event) 
  {
    /* reduce geometry and centroid bounds */
    size_t numAddedPrims = 0;
    BBox3f geomBound = empty;
    BBox3f centBound = empty;
    for (size_t i=0; i<numTasks; i++) {
      geomBound = merge(geomBound,geomBounds[i]);
      centBound = merge(centBound,centBounds[i]);
      numAddedPrims += work[i].numPrims;
    }
    new (&pinfo) PrimInfo(numAddedPrims,geomBound,centBound,pinfo);
    
    /* reduce heuristic and find best split */
    Heuristic heuristic; 
    Heuristic::reduce(heuristics,numTasks,heuristic); 
    heuristic.best(split);
  }
  
  /*! explicit template instantiations */
  INSTANTIATE_TEMPLATE_BY_HEURISTIC(PrimRefGen);
}
