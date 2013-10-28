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

#include "splitter_parallel.h"
#include "splitter_fallback.h"
#include "heuristics.h"

namespace embree
{
  template<typename Heuristic, typename PrimRefBlockList>
  MultiThreadedSplitter<Heuristic,PrimRefBlockList>::MultiThreadedSplitter(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event, 
                                                                           PrimRefAlloc* alloc, const RTCGeometry* geom,
                                                                           PrimRefBlockList& prims_i, const PrimInfo& pinfo, const Split& split, 
                                                                           TaskScheduler::completeFunction cfun, void* cptr)
    : alloc(alloc), prims(prims_i), pinfo(pinfo), split(split), geom(geom), cfun(cfun), cptr(cptr)
  {
    /* if split was not successfull enforce some split */
    if (unlikely(split.linfo.size() == 0 || split.rinfo.size() == 0)) {
      FallBackSplitter<Heuristic,PrimRefBlockList>::split(threadIndex,alloc,geom,prims,pinfo,lprims,linfo,lsplit,rprims,rinfo,rsplit);
      cfun(cptr,threadIndex,threadCount,event);
    } 
    /* perform spatial split */
    else if (unlikely(split.spatial())) 
    {
      new (&task) TaskScheduler::Task(event,
                                      _task_split_parallel_spatial,this,numTasks,
                                      _task_split_parallel_reduce,this,
                                      "build::parsplit");
      TaskScheduler::addTask(threadIndex,TaskScheduler::GLOBAL_FRONT,&task);
    }
    /* otherwise perform normal split */
    else {
      new (&task) TaskScheduler::Task(event,
                                      _task_split_parallel,       this,numTasks,
                                      _task_split_parallel_reduce,this,
                                      "build::parsplit");
      TaskScheduler::addTask(threadIndex,TaskScheduler::GLOBAL_FRONT,&task);
    }
  }

  template<typename Heuristic, typename PrimRefBlockList>
  void MultiThreadedSplitter<Heuristic,PrimRefBlockList>::task_split_parallel(size_t threadIndex, size_t threadCount, 
                                                                              size_t taskIndex, size_t taskCount,
                                                                              TaskScheduler::Event* event)
  {
    Heuristic& lheuristic = lheuristics[taskIndex]; new (&lheuristic) Heuristic (split.linfo,geom);
    Heuristic& rheuristic = rheuristics[taskIndex]; new (&rheuristic) Heuristic (split.rinfo,geom);
    atomic_set<PrimRefBlock>::item* lblock = lprims.insert(alloc->malloc(threadIndex));
    atomic_set<PrimRefBlock>::item* rblock = rprims.insert(alloc->malloc(threadIndex));
    
    while (atomic_set<PrimRefBlock>::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const PrimRef& prim = block->at(i); 

        if (split.left(prim)) 
        {
          if (likely(lblock->insert(prim))) continue; 
          lheuristic.bin(lblock->base(),lblock->size());
          lblock = lprims.insert(alloc->malloc(threadIndex));
          lblock->insert(prim);
        } 
        else 
        {
          if (likely(rblock->insert(prim))) continue;
          rheuristic.bin(rblock->base(),rblock->size());
          rblock = rprims.insert(alloc->malloc(threadIndex));
          rblock->insert(prim);
        }
      }
      alloc->free(threadIndex,block);
    }
    lheuristic.bin(lblock->base(),lblock->size());
    rheuristic.bin(rblock->base(),rblock->size());
  }

  template<typename Heuristic, typename PrimRefBlockList>
  void MultiThreadedSplitter<Heuristic,PrimRefBlockList>::task_split_parallel_spatial(size_t threadIndex, size_t threadCount, 
                                                                                      size_t taskIndex, size_t taskCount,
                                                                                      TaskScheduler::Event* event)
  {
    Heuristic& lheuristic = lheuristics[taskIndex]; new (&lheuristic) Heuristic (split.linfo,geom);
    Heuristic& rheuristic = rheuristics[taskIndex]; new (&rheuristic) Heuristic (split.rinfo,geom);
    atomic_set<PrimRefBlock>::item* lblock = lprims.insert(alloc->malloc(threadIndex));
    atomic_set<PrimRefBlock>::item* rblock = rprims.insert(alloc->malloc(threadIndex));
    
    while (atomic_set<PrimRefBlock>::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const PrimRef& prim = block->at(i); 
        PrimRef lprim, rprim; split.split(prim,lprim,rprim);

        if (lprim.id() != size_t(-1) && !lblock->insert(lprim)) 
        {
          lheuristic.bin(lblock->base(),lblock->size());
          lblock = lprims.insert(alloc->malloc(threadIndex));
          lblock->insert(lprim);
        } 
        
        if (rprim.id() != size_t(-1) && !rblock->insert(rprim)) 
        {
          rheuristic.bin(rblock->base(),rblock->size());
          rblock = rprims.insert(alloc->malloc(threadIndex));
          rblock->insert(rprim);
        }
      }
      alloc->free(threadIndex,block);
    }
    lheuristic.bin(lblock->base(),lblock->size());
    rheuristic.bin(rblock->base(),rblock->size());
  }

  template<typename Heuristic, typename PrimRefBlockList>
  void MultiThreadedSplitter<Heuristic,PrimRefBlockList>::task_split_parallel_reduce(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
  {
    Heuristic lheuristic; Heuristic::reduce(lheuristics,numTasks,lheuristic); linfo = split.linfo; lheuristic.best(lsplit); 
    Heuristic rheuristic; Heuristic::reduce(rheuristics,numTasks,rheuristic); rinfo = split.rinfo; rheuristic.best(rsplit); 
    assert(linfo.size());
    assert(rinfo.size());
    cfun(cptr,threadIndex,threadCount,event);
  }
  
  /*! explicit template instantiations */
  INSTANTIATE_TEMPLATE_BY_HEURISTIC_AND_PRIMREFLIST(MultiThreadedSplitter);
}
