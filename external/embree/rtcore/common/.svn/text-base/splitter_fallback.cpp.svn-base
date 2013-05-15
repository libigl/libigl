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

#include "splitter_fallback.h"
#include "heuristics.h"

namespace embree
{
  template<typename Heuristic, typename PrimRefBlockList>
  void FallBackSplitter<Heuristic,PrimRefBlockList>::split(const TaskScheduler::ThreadInfo& thread, PrimRefAlloc* alloc, 
                                                           const BuildTriangle* triangles, const Vec3fa* vertices,
                                                           PrimRefBlockList& prims, const PrimInfo& pinfo,
                                                           PrimRefBlockList& lprims, PrimInfo& linfo, Split& lsplit,
                                                           PrimRefBlockList& rprims, PrimInfo& rinfo, Split& rsplit)
  {
    /* enforce split */
    size_t lnum = 0; BBox3f lgeomBounds = empty; BBox3f lcentBounds = empty;
    size_t rnum = 0; BBox3f rgeomBounds = empty; BBox3f rcentBounds = empty;
    atomic_set<PrimRefBlock>::item* lblock = lprims.insert(alloc->malloc(thread));
    atomic_set<PrimRefBlock>::item* rblock = rprims.insert(alloc->malloc(thread));
    
    while (atomic_set<PrimRefBlock>::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const PrimRef& prim = block->at(i); 
        const BBox3f bounds = prim.bounds();
        
        if ((lnum+rnum)&1) 
        {
          lnum++;
          lgeomBounds.grow(bounds);
          lcentBounds.grow(center2(bounds));
          if (likely(lblock->insert(prim))) continue; 
          lblock = lprims.insert(alloc->malloc(thread));
          lblock->insert(prim);
        } else {
          rnum++;
          rgeomBounds.grow(bounds);
          rcentBounds.grow(center2(bounds));
          if (likely(rblock->insert(prim))) continue;
          rblock = rprims.insert(alloc->malloc(thread));
          rblock->insert(prim);
        }
      }
    }
    new (&linfo) PrimInfo(lnum,lgeomBounds,lcentBounds);
    new (&rinfo) PrimInfo(rnum,rgeomBounds,rcentBounds);

    /* perform binning of left side */
    Heuristic lheuristic(linfo,triangles,vertices);
    typename PrimRefBlockList::iterator liter(lprims);
    while (typename PrimRefBlockList::item* block = liter.next()) {
      lheuristic.bin(block->base(),block->size());
    }
    lheuristic.best(lsplit);

    /* perform binning of right side */
    Heuristic rheuristic(rinfo,triangles,vertices);
    typename PrimRefBlockList::iterator riter(rprims);
    while (typename PrimRefBlockList::item* block = riter.next()) {
      rheuristic.bin(block->base(),block->size());
    }
    rheuristic.best(rsplit);
  }

  /*! explicit template instantiations */
  INSTANTIATE_TEMPLATE_BY_HEURISTIC_AND_PRIMREFLIST(FallBackSplitter);
}
