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

#include "splitter_parallel.h"
#include "splitter_fallback.h"
#include "heuristics.h"

namespace embree
{
  template<typename Heuristic, typename PrimRefBlockList>
  MultiThreadedSplitter<Heuristic,PrimRefBlockList>::MultiThreadedSplitter(const TaskScheduler::ThreadInfo& thread,
                                                                           PrimRefAlloc* alloc, const BuildTriangle* triangles, const Vec3fa* vertices,
                                                                           PrimRefBlockList& prims_i, const PrimInfo& pinfo, const Split& split, 
                                                                           TaskScheduler::completeFunction cfun, void* cptr)
    : alloc(alloc), prims(prims_i), pinfo(pinfo), split(split), triangles(triangles), vertices(vertices), cfun(cfun), cptr(cptr)
  {
    /* if split was not successfull enforce some split */
    if (unlikely(split.linfo.size() == 0 || split.rinfo.size() == 0)) {
      FallBackSplitter<Heuristic,PrimRefBlockList>::split(thread,alloc,triangles,vertices,prims,pinfo,lprims,linfo,lsplit,rprims,rinfo,rsplit);
      cfun(thread,cptr);
    } 
    /* perform spatial split */
    else if (unlikely(split.spatial())) 
    {
      scheduler->addTask(thread,TaskScheduler::GLOBAL_FRONT,
                         (TaskScheduler::runFunction     )_task_split_parallel_spatial,this,numTasks,
                         (TaskScheduler::completeFunction)_task_split_parallel_reduce,this,
                         "build::parsplit");
    }
    /* otherwise perform normal split */
    else {
      scheduler->addTask(thread,TaskScheduler::GLOBAL_FRONT,
                         (TaskScheduler::runFunction     )_task_split_parallel,       this,numTasks,
                         (TaskScheduler::completeFunction)_task_split_parallel_reduce,this,
                         "build::parsplit");
    }
  }

  template<typename Heuristic, typename PrimRefBlockList>
  void MultiThreadedSplitter<Heuristic,PrimRefBlockList>::task_split_parallel(const TaskScheduler::ThreadInfo& thread, size_t idx)
  {
    Heuristic& lheuristic = lheuristics[idx]; new (&lheuristic) Heuristic (split.linfo,triangles,vertices);
    Heuristic& rheuristic = rheuristics[idx]; new (&rheuristic) Heuristic (split.rinfo,triangles,vertices);
    atomic_set<PrimRefBlock>::item* lblock = lprims.insert(alloc->malloc(thread));
    atomic_set<PrimRefBlock>::item* rblock = rprims.insert(alloc->malloc(thread));
    
    while (atomic_set<PrimRefBlock>::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const PrimRef& prim = block->at(i); 

        if (split.left(prim)) 
        {
          if (likely(lblock->insert(prim))) continue; 
          lheuristic.bin(lblock->base(),lblock->size());
          lblock = lprims.insert(alloc->malloc(thread));
          lblock->insert(prim);
        } 
        else 
        {
          if (likely(rblock->insert(prim))) continue;
          rheuristic.bin(rblock->base(),rblock->size());
          rblock = rprims.insert(alloc->malloc(thread));
          rblock->insert(prim);
        }
      }
      alloc->free(thread,block);
    }
    lheuristic.bin(lblock->base(),lblock->size());
    rheuristic.bin(rblock->base(),rblock->size());
  }

  template<typename Heuristic, typename PrimRefBlockList>
  void MultiThreadedSplitter<Heuristic,PrimRefBlockList>::task_split_parallel_spatial(const TaskScheduler::ThreadInfo& thread, size_t idx)
  {
    Heuristic& lheuristic = lheuristics[idx]; new (&lheuristic) Heuristic (split.linfo,triangles,vertices);
    Heuristic& rheuristic = rheuristics[idx]; new (&rheuristic) Heuristic (split.rinfo,triangles,vertices);
    atomic_set<PrimRefBlock>::item* lblock = lprims.insert(alloc->malloc(thread));
    atomic_set<PrimRefBlock>::item* rblock = rprims.insert(alloc->malloc(thread));
    
    while (atomic_set<PrimRefBlock>::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const PrimRef& prim = block->at(i); 
        PrimRef lprim, rprim; split.split(prim,lprim,rprim);

        if (lprim.id() != size_t(-1) && !lblock->insert(lprim)) 
        {
          lheuristic.bin(lblock->base(),lblock->size());
          lblock = lprims.insert(alloc->malloc(thread));
          lblock->insert(lprim);
        } 
        
        if (rprim.id() != size_t(-1) && !rblock->insert(rprim)) 
        {
          rheuristic.bin(rblock->base(),rblock->size());
          rblock = rprims.insert(alloc->malloc(thread));
          rblock->insert(rprim);
        }
      }
      alloc->free(thread,block);
    }
    lheuristic.bin(lblock->base(),lblock->size());
    rheuristic.bin(rblock->base(),rblock->size());
  }

  template<typename Heuristic, typename PrimRefBlockList>
  void MultiThreadedSplitter<Heuristic,PrimRefBlockList>::task_split_parallel_reduce(const TaskScheduler::ThreadInfo& thread)
  {
    Heuristic lheuristic; Heuristic::reduce(lheuristics,numTasks,lheuristic); linfo = split.linfo; lheuristic.best(lsplit); 
    Heuristic rheuristic; Heuristic::reduce(rheuristics,numTasks,rheuristic); rinfo = split.rinfo; rheuristic.best(rsplit); 
    assert(linfo.size());
    assert(rinfo.size());
    cfun(thread,cptr);
  }
  
  /*! explicit template instantiations */
  INSTANTIATE_TEMPLATE_BY_HEURISTIC_AND_PRIMREFLIST(MultiThreadedSplitter);
}
