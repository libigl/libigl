// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
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

#include "heuristic_strand_partition.h"

namespace embree
{
  namespace isa
  {
    template<>
    const StrandSplit::Split StrandSplit::find<false>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, BezierRefList& prims)
    {
      /* first curve determines first axis */
      BezierRefList::block_iterator_unsafe i = prims;
      Vec3fa axis0 = normalize(i->p3 - i->p0);
      
      /* find 2nd axis that is most misaligned with first axis */
      float bestCos = 1.0f;
      Vec3fa axis1 = axis0;
      for (; i; i++) {
	Vec3fa axisi = i->p3 - i->p0;
	float leni = length(axisi);
	if (leni == 0.0f) continue;
	axisi /= leni;
	float cos = abs(dot(axisi,axis0));
	if (cos < bestCos) { bestCos = cos; axis1 = axisi; }
      }
      
      /* partition the two strands */
      size_t lnum = 0, rnum = 0;
      BBox3fa lbounds = empty, rbounds = empty;
      const LinearSpace3fa space0 = frame(axis0).transposed();
      const LinearSpace3fa space1 = frame(axis1).transposed();
      
      for (BezierRefList::block_iterator_unsafe i = prims; i; i++) 
      {
	BezierPrim& prim = *i;
	const Vec3fa axisi = normalize(prim.p3-prim.p0);
	const float cos0 = abs(dot(axisi,axis0));
	const float cos1 = abs(dot(axisi,axis1));
	
	if (cos0 > cos1) { lnum++; lbounds.extend(prim.bounds(space0)); }
	else             { rnum++; rbounds.extend(prim.bounds(space1)); }
      }
      
      /*! return an invalid split if we do not partition */
      if (lnum == 0 || rnum == 0) 
	return Split(inf,axis0,axis1);
      
      /*! calculate sah for the split */
      const float sah = float(lnum)*halfArea(lbounds) + float(rnum)*halfArea(rbounds);
      return Split(sah,axis0,axis1);
    }
    
    StrandSplit::TaskFindParallel::TaskFindParallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, BezierRefList& prims)
    {
      /* first curve determines first axis */
      BezierRefList::block_iterator_unsafe i = prims;
      axis0 = axis1 = normalize(i->p3 - i->p0);
      
      /* parallel calculation of 2nd axis  */
      size_t numTasks = min(maxTasks,threadCount);
      scheduler->dispatchTask(threadIndex,numTasks,_task_bound_parallel,this,numTasks,"build::task_find_parallel");
      
      /* select best 2nd axis */
      float bestCos = 1.0f;
      for (size_t i=0; i<numTasks; i++) {
	if (task_cos[i] < bestCos) { bestCos = task_cos[i]; axis1 = task_axis1[i]; }
      }
      
      /* parallel calculation of unaligned bounds */
      scheduler->dispatchTask(threadIndex,numTasks,_task_bound_parallel,this,numTasks,"build::task_find_parallel");
      
      /* reduce bounds calculates by tasks */
      size_t lnum = 0; BBox3fa lbounds = empty;
      size_t rnum = 0; BBox3fa rbounds = empty;
      for (size_t i=0; i<numTasks; i++) {
	lnum += task_lnum[i]; lbounds.extend(task_lbounds[i]);
	rnum += task_rnum[i]; rbounds.extend(task_rbounds[i]);
      }
      
      /*! return an invalid split if we do not partition */
      if (lnum == 0 || rnum == 0) {
	split = Split(inf,axis0,axis1); 
	return;
      }
      
      /*! calculate sah for the split */
      const float sah = float(lnum)*halfArea(lbounds) + float(lnum)*halfArea(rbounds);
      split = Split(sah,axis0,axis1);
    }
    
    void StrandSplit::TaskFindParallel::task_find_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      float bestCos = 1.0f;
      Vec3fa bestAxis1 = axis0;
      
      while (BezierRefList::item* block = iter0.next()) 
      {
	for (size_t i=0; i<block->size(); i++) 
	{
	  BezierPrim& prim = block->at(i);
	  Vec3fa axisi = prim.p3 - prim.p0;
	  float leni = length(axisi);
	  if (leni == 0.0f) continue;
	  axisi /= leni;
	  float cos = abs(dot(axisi,axis0));
	  if (cos < bestCos) { bestCos = cos; bestAxis1 = axisi; }
	}
      }
      task_cos[taskIndex] = bestCos;
      task_axis1[taskIndex] = bestAxis1;
    }
    
    void StrandSplit::TaskFindParallel::task_bound_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      const LinearSpace3fa space0 = frame(axis0).transposed();
      const LinearSpace3fa space1 = frame(axis1).transposed();
      
      size_t lnum = 0; BBox3fa lbounds = empty;
      size_t rnum = 0; BBox3fa rbounds = empty;
      while (BezierRefList::item* block = iter1.next()) 
      {
	for (size_t i=0; i<block->size(); i++) 
	{
	  BezierPrim& prim = block->at(i);
	  const Vec3fa axisi = normalize(prim.p3-prim.p0);
	  const float cos0 = abs(dot(axisi,axis0));
	  const float cos1 = abs(dot(axisi,axis1));
	  
	  if (cos0 > cos1) { lnum++; lbounds.extend(prim.bounds(space0)); }
	  else             { rnum++; rbounds.extend(prim.bounds(space1)); }
	}
      }
      task_lnum[taskIndex] = lnum; task_lbounds[taskIndex] = lbounds;
      task_rnum[taskIndex] = rnum; task_rbounds[taskIndex] = rbounds;
    }
    
    template<>
    const StrandSplit::Split StrandSplit::find<true>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, BezierRefList& prims) {
      return TaskFindParallel(threadIndex,threadCount,scheduler,prims).split;
    }
    
    template<>
    void StrandSplit::Split::split<false>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefBlockAlloc<BezierPrim>& alloc, 
					  BezierRefList& prims, 
					  BezierRefList& lprims_o, PrimInfo& linfo_o, 
					  BezierRefList& rprims_o, PrimInfo& rinfo_o) const 
    {
      BezierRefList::item* lblock = lprims_o.insert(alloc.malloc(threadIndex));
      BezierRefList::item* rblock = rprims_o.insert(alloc.malloc(threadIndex));
      linfo_o.reset();
      rinfo_o.reset();
      
      while (BezierRefList::item* block = prims.take()) 
      {
	for (size_t i=0; i<block->size(); i++) 
	{
	  const BezierPrim& prim = block->at(i); 
	  const Vec3fa axisi = normalize(prim.p3-prim.p0);
	  const float cos0 = abs(dot(axisi,axis0));
	  const float cos1 = abs(dot(axisi,axis1));
	  
	  if (cos0 > cos1) 
	  {
	    linfo_o.add(prim.bounds(),prim.center());
	    if (likely(lblock->insert(prim))) continue; 
	    lblock = lprims_o.insert(alloc.malloc(threadIndex));
	    lblock->insert(prim);
	  } 
	  else 
	  {
	    rinfo_o.add(prim.bounds(),prim.center());
	    if (likely(rblock->insert(prim))) continue;
	    rblock = rprims_o.insert(alloc.malloc(threadIndex));
	    rblock->insert(prim);
	  }
	}
	alloc.free(threadIndex,block);
      }
    }
    
    StrandSplit::TaskSplitParallel::TaskSplitParallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, const Split* split, PrimRefBlockAlloc<BezierPrim>& alloc, 
						      BezierRefList& prims, 
						      BezierRefList& lprims_o, PrimInfo& linfo_o, 
						      BezierRefList& rprims_o, PrimInfo& rinfo_o)
      : split(split), alloc(alloc), prims(prims), lprims_o(lprims_o), linfo_o(linfo_o), rprims_o(rprims_o), rinfo_o(rinfo_o)
    {
      /* parallel calculation of centroid bounds */
      size_t numTasks = min(maxTasks,threadCount);
      scheduler->dispatchTask(threadIndex,numTasks,_task_split_parallel,this,numTasks,"build::task_split_parallel");
      
      /* reduction of bounding info */
      linfo_o = linfos[0];
      rinfo_o = rinfos[0];
      for (size_t i=1; i<numTasks; i++) {
	linfo_o.merge(linfos[i]);
	rinfo_o.merge(rinfos[i]);
      }
    }
    
    void StrandSplit::TaskSplitParallel::task_split_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      split->split<false>(threadIndex,threadCount,NULL,alloc,prims,lprims_o,linfos[taskIndex],rprims_o,rinfos[taskIndex]);
    }
    
    template<>
    void StrandSplit::Split::split<true>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, 
					 PrimRefBlockAlloc<BezierPrim>& alloc, BezierRefList& prims, 
					 BezierRefList& lprims_o, PrimInfo& linfo_o, 
					 BezierRefList& rprims_o, PrimInfo& rinfo_o) const
    {
      TaskSplitParallel(threadIndex,threadCount,scheduler,this,alloc,prims,lprims_o,linfo_o,rprims_o,rinfo_o);
    }
  }
}
