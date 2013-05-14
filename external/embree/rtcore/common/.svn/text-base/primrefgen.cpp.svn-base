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

#include "primrefgen.h"
#include "heuristics.h"

namespace embree
{
  template<typename Heuristic, typename PrimRefBlockList>
  PrimRefGen<Heuristic,PrimRefBlockList>::PrimRefGen(const TaskScheduler::ThreadInfo& thread,
                                                     const BuildTriangle* triangles, size_t numTriangles, 
                                                     const Vec3fa* vertices, size_t numVertices,
                                                     const BBox3f& bounds, 
                                                     PrimRefAlloc* alloc)
    
    : triangles(triangles), numTriangles(numTriangles), vertices(vertices), numVertices(numVertices), alloc(alloc)
  {
    /* approximate bounds if not available */
    if (bounds.empty() && numTriangles) {
      BBox3f geomBounds = empty, centBounds = empty;
      for (size_t i=0; i<2048; i++) {
        MakePrimRef(i*numTriangles/2048,geomBounds,centBounds);
      }
      new (&pinfo) PrimInfo(numTriangles,geomBounds,centBounds);
    }
    else
      new (&pinfo) PrimInfo(numTriangles,bounds);
    
    /* start parallel task */
    scheduler->addTask(thread,TaskScheduler::GLOBAL_FRONT,
                       (TaskScheduler::runFunction     )_task_gen_parallel,       this,numTasks,
                       (TaskScheduler::completeFunction)_task_gen_parallel_reduce,this,
                       "build::primrefgen");
  }

  template<typename Heuristic, typename PrimRefBlockList>
  PrimRefGen<Heuristic,PrimRefBlockList>::~PrimRefGen() {
    pinfo.clear();
  }

  template<typename Heuristic, typename PrimRefBlockList>
  __forceinline PrimRef PrimRefGen<Heuristic,PrimRefBlockList>::MakePrimRef(size_t i, BBox3f& geomBound, BBox3f& centBound) 
  {
    const BuildTriangle& tri = triangles[i];
    const Vec3f& v0 = vertices[tri.v0];
    const Vec3f& v1 = vertices[tri.v1];
    const Vec3f& v2 = vertices[tri.v2];
    const BBox3f b = merge(BBox3f(v0), BBox3f(v1), BBox3f(v2));
    geomBound.grow(b);
    centBound.grow(center2(b));
    return PrimRef(b,i);
  }

  template<typename Heuristic, typename PrimRefBlockList>
  void PrimRefGen<Heuristic,PrimRefBlockList>::task_gen_parallel(const TaskScheduler::ThreadInfo& thread, size_t idx) 
  {
    Heuristic& heuristic = heuristics[idx];
    new (&heuristic) Heuristic(pinfo,triangles,vertices);
    
    /* static work allocation */
    size_t start = (idx+0)*numTriangles/numTasks;
    size_t end   = (idx+1)*numTriangles/numTasks;
    BBox3f geomBound = empty, centBound = empty;
    typename PrimRefBlockList::item* block = prims.insert(alloc->malloc(thread)); 
    
    for (size_t i=start; i<end; i++) 
    {
      const PrimRef prim = MakePrimRef(i,geomBound,centBound);
      if (likely(block->insert(prim))) continue; 
      heuristic.bin(block->base(),block->size());
      block = prims.insert(alloc->malloc(thread));
      block->insert(prim);
    }
    heuristic.bin(block->base(),block->size());
    geomBounds[idx] = geomBound;
    centBounds[idx] = centBound;
  }
  
  template<typename Heuristic, typename PrimRefBlockList>
  void PrimRefGen<Heuristic,PrimRefBlockList>::task_gen_parallel_reduce(const TaskScheduler::ThreadInfo& thread) 
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
