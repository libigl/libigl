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

#include "splitter.h"
#include "splitter_fallback.h"
#include "heuristics.h"

namespace embree
{
  template<typename Heuristic>
  Splitter<Heuristic>::Splitter (const TaskScheduler::ThreadInfo& thread, PrimRefAlloc* alloc, const BuildTriangle* triangles, const Vec3fa* vertices,
                                 atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& psplit)
  {
    /* if split was not successfull enforce some split */
    if (unlikely(psplit.linfo.size() == 0 || psplit.rinfo.size() == 0)) {
      FallBackSplitter<Heuristic,atomic_set<PrimRefBlock> >::split(thread,alloc,triangles,vertices,prims,pinfo,lprims,linfo,lsplit,rprims,rinfo,rsplit);
    }
    /* split with support for duplications */
    else if (unlikely(psplit.spatial())) {
      split_spatial(thread,alloc,triangles,vertices,prims,pinfo,psplit);
    }
    /* otherwise perform normal split */
    else {
      split(thread,alloc,triangles,vertices,prims,pinfo,psplit);
    }
    assert(linfo.size());
    assert(rinfo.size());
  }

  template<typename Heuristic>
  void Splitter<Heuristic>::split(const TaskScheduler::ThreadInfo& thread, PrimRefAlloc* alloc, const BuildTriangle* triangles, const Vec3fa* vertices,
                                  atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
  {
    Heuristic lheuristic(split.linfo,triangles,vertices);
    Heuristic rheuristic(split.rinfo,triangles,vertices);
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
    lheuristic.bin(lblock->base(),lblock->size()); linfo = split.linfo; lheuristic.best(lsplit); 
    rheuristic.bin(rblock->base(),rblock->size()); rinfo = split.rinfo; rheuristic.best(rsplit);
  }
  
  template<typename Heuristic>
  void Splitter<Heuristic>::split_spatial(const TaskScheduler::ThreadInfo& thread, PrimRefAlloc* alloc, const BuildTriangle* triangles, const Vec3fa* vertices,
                                          atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
  {
    Heuristic lheuristic(split.linfo,triangles,vertices);
    Heuristic rheuristic(split.rinfo,triangles,vertices);
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
    lheuristic.bin(lblock->base(),lblock->size()); linfo = split.linfo; lheuristic.best(lsplit); 
    rheuristic.bin(rblock->base(),rblock->size()); rinfo = split.rinfo; rheuristic.best(rsplit);
  }

  /*! explicit template instantiations */
  INSTANTIATE_TEMPLATE_BY_HEURISTIC(Splitter);
}
