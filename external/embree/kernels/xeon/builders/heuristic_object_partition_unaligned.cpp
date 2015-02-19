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

#include "heuristic_object_partition_unaligned.h"

namespace embree
{
  namespace isa
  {
    const LinearSpace3fa ObjectPartitionUnaligned::computeAlignedSpace(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, BezierRefList& prims)
    {
      /*! find first curve that defines valid direction */
      Vec3fa axis(0,0,1);
      BezierRefList::block_iterator_unsafe i = prims;
      for (; i; i++)
      {
	const Vec3fa axis1 = normalize(i->p3 - i->p0);
	if (length(i->p3 - i->p0) > 1E-9f) {
	  axis = axis1;
	  break;
	}
      }
      return frame(axis).transposed();
    }
    
    const std::pair<AffineSpace3fa,AffineSpace3fa> ObjectPartitionUnaligned::computeAlignedSpaceMB(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, Scene* scene, BezierRefList& prims)
    {
      /*! find first curve that defines valid direction */
      Vec3fa p0(0,0,0);
      Vec3fa p1(0,0,0);
      Vec3fa axis0(0,0,1), axisb0(0,1,0);
      Vec3fa axis1(0,0,1), axisb1(0,1,0);
      //Vec3fa axis2(0,1,0);
      BezierRefList::block_iterator_unsafe i = prims;
      for (; i; i++)
      {
        const BezierCurves* curves = scene->getBezierCurves(i->geomID<0>());
        const int curve = curves->curve(i->primID<0>());

	const Vec3fa a3 = curves->vertex(curve+3,0);
	const Vec3fa a2 = curves->vertex(curve+2,0);
	const Vec3fa a1 = curves->vertex(curve+1,0);
        const Vec3fa a0 = curves->vertex(curve+0,0);
					 
        const Vec3fa b3 = curves->vertex(curve+3,1);
	const Vec3fa b2 = curves->vertex(curve+2,1);
	const Vec3fa b1 = curves->vertex(curve+1,1);
        const Vec3fa b0 = curves->vertex(curve+0,1);

	if (length(a3 - a0) > 1E-9f && length(a1 - a0) > 1E-9f &&
	    length(b3 - b0) > 1E-9f && length(b1 - b0) > 1E-9f) 
	{
	  axis0 = normalize(a3 - a0); axisb0 = normalize(a1 - a0); 
	  axis1 = normalize(b3 - b0); axisb1 = normalize(b1 - b0); 
	  p0 = a0; p1 = b0;
          /*if (length(b3-a3) > 1E-9f) axis2 = b3-a3;
          else if (length(b0-a0) > 1E-9f) axis2 = b0-a0;
          else axis2 = Vec3fa(1,0,0);*/ // FIXME: not correct
	  break;
	}
      }
#if 1
      //LinearSpace3fa space01 = frame(0.5f*axis0 + 0.5f*axis1).transposed();
      //AffineSpace3fa space0 = frame(axis0).transposed();
      //AffineSpace3fa space1 = frame(axis1).transposed();

      const Vec3fa space0_dx = normalize(axis0);
      const Vec3fa space0_dy = normalize(cross(space0_dx,axisb0));
      const Vec3fa space0_dz = normalize(cross(space0_dx,space0_dy));
      LinearSpace3fa space0(space0_dz,space0_dy,space0_dx);
      space0 = space0.transposed();

      const Vec3fa space1_dx = normalize(axis1);
      const Vec3fa space1_dy = normalize(cross(space1_dx,axisb1));
      const Vec3fa space1_dz = normalize(cross(space1_dx,space1_dy));
      LinearSpace3fa space1(space1_dz,space1_dy,space1_dx);
      space1 = space1.transposed();

      //space0.p = -xfmVector(space0.l,p0);
      //space1.p = -xfmVector(space1.l,p1);
      //space0 = space01;
      //space1 = space0;
#else
      const Vec3fa space0_dx = normalize(axis0);
      const Vec3fa space0_dy = normalize(cross(space0_dx,axis2));
      const Vec3fa space0_dz = normalize(cross(space0_dx,space0_dy));
      AffineSpace3fa space0(space0_dz,space0_dy,space0_dx,zero);
      space0.l = space0.l.transposed();

      const Vec3fa space1_dx = normalize(axis1);
      const Vec3fa space1_dy = normalize(cross(space1_dx,axis2));
      const Vec3fa space1_dz = normalize(cross(space1_dx,space1_dy));
      AffineSpace3fa space1(space1_dz,space1_dy,space1_dx,zero);
      space1.l = space1.l.transposed();

      space0.p = -xfmVector(space0.l,p0);
      space1.p = -xfmVector(space1.l,p1);
#endif 

      return std::pair<AffineSpace3fa,AffineSpace3fa>(space0,space1);
    }

    template<>
    const PrimInfo ObjectPartitionUnaligned::computePrimInfo<false>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, BezierRefList& prims, const LinearSpace3fa& space)
    {
      size_t num = 0;
      BBox3fa geomBounds = empty;
      BBox3fa centBounds = empty;
      for (BezierRefList::block_iterator_unsafe i = prims; i; i++) {
	num++;
	geomBounds.extend(i->bounds(space));
	centBounds.extend(i->center(space));
      }
      return PrimInfo(num,geomBounds,centBounds);
    }
    
    template<>
    const PrimInfo ObjectPartitionUnaligned::computePrimInfo<true>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, BezierRefList& prims, const LinearSpace3fa& space)
    {
      const TaskPrimInfoParallel bounds(threadIndex,threadCount,scheduler,prims,space);
      return PrimInfo(bounds.num,bounds.geomBounds,bounds.centBounds);
    }

    template<>
    const ObjectPartitionUnaligned::PrimInfoMB ObjectPartitionUnaligned::computePrimInfoMB<false>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, Scene* scene, BezierRefList& prims, 
                                                                        const std::pair<AffineSpace3fa,AffineSpace3fa>& spaces)
    {
      size_t N = 0;
      BBox3fa centBounds = empty;
      BBox3fa geomBounds = empty;
      BBox3fa s0t0 = empty, s0t1_s1t0 = empty, s1t1 = empty;
      for (BezierRefList::block_iterator_unsafe i = prims; i; i++) 
      {
	N++;
        geomBounds.extend(i->bounds(spaces.first));
        centBounds.extend(i->center(spaces.first));

        const BezierCurves* curves = scene->getBezierCurves(i->geomID<0>());
	s0t0.extend(curves->bounds(spaces.first,i->primID<0>(),0));
        s0t1_s1t0.extend(curves->bounds(spaces.first,spaces.second,i->primID<0>()));
	s1t1.extend(curves->bounds(spaces.second,i->primID<0>(),1));
      }

      PrimInfoMB ret;
      ret.pinfo = PrimInfo(N,geomBounds,centBounds);
      ret.s0t0 = s0t0;
      ret.s0t1_s1t0 = s0t1_s1t0;
      ret.s1t1 = s1t1;
      return ret;
    }
    
    template<>
    const ObjectPartitionUnaligned::PrimInfoMB ObjectPartitionUnaligned::computePrimInfoMB<true>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, Scene* scene, BezierRefList& prims, 
                                                                       const std::pair<AffineSpace3fa,AffineSpace3fa>& spaces)
    {
      const TaskPrimInfoMBParallel bounds(threadIndex,threadCount,scheduler,scene,prims,spaces.first,spaces.second);

      PrimInfoMB ret;
      ret.pinfo = PrimInfo(bounds.num,bounds.geomBounds,bounds.centBounds);
      ret.s0t0 = bounds.s0t0;
      ret.s0t1_s1t0 = bounds.s0t1_s1t0;
      ret.s1t1 = bounds.s1t1;
      return ret;
    }
    
    __forceinline ObjectPartitionUnaligned::Mapping::Mapping(const BBox3fa& centBounds, const LinearSpace3fa& space) 
      : space(space)
    {
      const ssef diag = (ssef) centBounds.size();
      scale = select(diag > ssef(1E-19f),rcp(diag) * ssef(BINS * 0.99f),ssef(0.0f));
      ofs  = (ssef) centBounds.lower;
    }
    
    __forceinline ssei ObjectPartitionUnaligned::Mapping::bin(const Vec3fa& p) const 
    {
      const ssei i = floori((ssef(p)-ofs)*scale);
#if 1
      assert(i[0] >=0 && i[0] < BINS);
      assert(i[1] >=0 && i[1] < BINS);
      assert(i[2] >=0 && i[2] < BINS);
      return i;
#else
      return clamp(i,ssei(0),ssei(BINS-1));
#endif
    }
    
    __forceinline ssei ObjectPartitionUnaligned::Mapping::bin_unsafe(const Vec3fa& p) const {
      return floori((ssef(p)-ofs)*scale);
    }
    
    __forceinline bool ObjectPartitionUnaligned::Mapping::invalid(const int dim) const {
      return scale[dim] == 0.0f;
    }
    
    __forceinline ObjectPartitionUnaligned::BinInfo::BinInfo() 
    {
      for (size_t i=0; i<BINS; i++) {
	bounds[i][0] = bounds[i][1] = bounds[i][2] = bounds[i][3] = empty;
	counts[i] = 0;
      }
    }
    
    __forceinline void ObjectPartitionUnaligned::BinInfo::bin (const BezierPrim* prims, size_t N, const Mapping& mapping)
    {
      for (size_t i=0; i<N; i++)
      {
	const BBox3fa cbounds = prims[i].bounds(mapping.space);
	const Vec3fa  center  = prims[i].center(mapping.space);
	const ssei bin = mapping.bin(center);
	const int b0 = bin[0]; counts[b0][0]++; bounds[b0][0].extend(cbounds);
	const int b1 = bin[1]; counts[b1][1]++; bounds[b1][1].extend(cbounds);
	const int b2 = bin[2]; counts[b2][2]++; bounds[b2][2].extend(cbounds);
      }
    }
    
    __forceinline void  ObjectPartitionUnaligned::BinInfo::bin(BezierRefList& prims, const Mapping& mapping)
    {
      BezierRefList::iterator i=prims;
      while (BezierRefList::item* block = i.next())
	bin(block->base(),block->size(),mapping);
    }
    
    __forceinline void ObjectPartitionUnaligned::BinInfo::merge (const BinInfo& other)
    {
      for (size_t i=0; i<BINS; i++)  // FIXME: dont iterate over all bins
      {
	counts[i] += other.counts[i];
	bounds[i][0].extend(other.bounds[i][0]);
	bounds[i][1].extend(other.bounds[i][1]);
	bounds[i][2].extend(other.bounds[i][2]);
      }
    }
    
    __forceinline ObjectPartitionUnaligned::Split ObjectPartitionUnaligned::BinInfo::best(BezierRefList& prims, const Mapping& mapping)
    {
      /* sweep from right to left and compute parallel prefix of merged bounds */
      ssef rAreas[BINS];
      ssei rCounts[BINS];
      ssei count = 0; BBox3fa bx = empty; BBox3fa by = empty; BBox3fa bz = empty;
      for (size_t i=BINS-1; i>0; i--)
      {
	count += counts[i];
	rCounts[i] = count;
	bx.extend(bounds[i][0]); rAreas[i][0] = halfArea(bx);
	by.extend(bounds[i][1]); rAreas[i][1] = halfArea(by);
	bz.extend(bounds[i][2]); rAreas[i][2] = halfArea(bz);
      }
      
      /* sweep from left to right and compute SAH */
      ssei ii = 1; ssef vbestSAH = pos_inf; ssei vbestPos = 0;
      count = 0; bx = empty; by = empty; bz = empty;
      for (size_t i=1; i<BINS; i++, ii+=1)
      {
	count += counts[i-1];
	bx.extend(bounds[i-1][0]); float Ax = halfArea(bx);
	by.extend(bounds[i-1][1]); float Ay = halfArea(by);
	bz.extend(bounds[i-1][2]); float Az = halfArea(bz);
	const ssef lArea = ssef(Ax,Ay,Az,Az);
	const ssef rArea = rAreas[i];
	const ssei lCount = blocks(count);
	const ssei rCount = blocks(rCounts[i]);
	const ssef sah = lArea*ssef(lCount) + rArea*ssef(rCount);
	vbestPos = select(sah < vbestSAH,ii ,vbestPos);
	vbestSAH = select(sah < vbestSAH,sah,vbestSAH);
      }
      
      /* find best dimension */
      float bestSAH = inf;
      int   bestDim = -1;
      int   bestPos = 0;
      for (size_t dim=0; dim<3; dim++) 
      {
	/* ignore zero sized dimensions */
	if (unlikely(mapping.invalid(dim)))
	  continue;
	
	/* test if this is a better dimension */
	if (vbestSAH[dim] < bestSAH && vbestPos[dim] != 0) {
	  bestDim = dim;
	  bestPos = vbestPos[dim];
	  bestSAH = vbestSAH[dim];
	}
      }
      
      return ObjectPartitionUnaligned::Split(bestSAH,bestDim,bestPos,mapping);
    }
    
    template<>
    const ObjectPartitionUnaligned::Split ObjectPartitionUnaligned::find<false>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, BezierRefList& prims, const LinearSpace3fa& space, const PrimInfo& pinfo)
    {
      BinInfo binner;
      const Mapping mapping(pinfo.centBounds,space);
      binner.bin(prims,mapping);
      return binner.best(prims,mapping);
    }

    ObjectPartitionUnaligned::TaskPrimInfoMBParallel::TaskPrimInfoMBParallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, Scene* scene, BezierRefList& prims, const AffineSpace3fa& space0, const AffineSpace3fa& space1) 
      : scene(scene), space0(space0), space1(space1), iter(prims), geomBounds(empty), centBounds(empty), s0t0(empty), s0t1_s1t0(empty), s1t1(empty)
    {
      size_t numTasks = min(maxTasks,threadCount);
      scheduler->dispatchTask(threadIndex,numTasks,_task_bound_parallel,this,numTasks,"build::task_bound_parallel");
    }
    
    void ObjectPartitionUnaligned::TaskPrimInfoMBParallel::task_bound_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      size_t N = 0;
      BBox3fa centBounds = empty;
      BBox3fa geomBounds = empty;
      BBox3fa s0t0 = empty, s0t1_s1t0 = empty, s1t1 = empty;
      while (BezierRefList::item* block = iter.next()) 
      {
	N += block->size();
	for (size_t i=0; i<block->size(); i++) 
        {
	  geomBounds.extend(block->at(i).bounds(space0));
	  centBounds.extend(block->at(i).center(space0));

          const BezierPrim& ref = block->at(i);
          const BezierCurves* curves = scene->getBezierCurves(ref.geomID<0>());
          s0t0.extend(curves->bounds(space0,ref.primID<0>(),0));
          s0t1_s1t0.extend(curves->bounds(space0,space1,ref.primID<0>()));
	  s1t1.extend(curves->bounds(space1,ref.primID<0>(),1));
	}
      }
      atomic_add(&this->num,N);
      this->centBounds.extend_atomic(centBounds);
      this->geomBounds.extend_atomic(geomBounds);
      this->s0t0.extend_atomic(s0t0);
      this->s0t1_s1t0.extend_atomic(s0t1_s1t0);
      this->s1t1.extend_atomic(s1t1);
    }
    
    ObjectPartitionUnaligned::TaskPrimInfoParallel::TaskPrimInfoParallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, BezierRefList& prims, const LinearSpace3fa& space) 
      : space(space), iter(prims), geomBounds(empty), centBounds(empty)
    {
      size_t numTasks = min(maxTasks,threadCount);
      scheduler->dispatchTask(threadIndex,numTasks,_task_bound_parallel,this,numTasks,"build::task_bound_parallel");
    }
    
    void ObjectPartitionUnaligned::TaskPrimInfoParallel::task_bound_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      size_t N = 0;
      BBox3fa centBounds = empty;
      BBox3fa geomBounds = empty;
      while (BezierRefList::item* block = iter.next()) 
      {
	N += block->size();
	for (size_t i=0; i<block->size(); i++) {
	  geomBounds.extend(block->at(i).bounds(space));
	  centBounds.extend(block->at(i).center(space));
	}
      }
      atomic_add(&this->num,N);
      this->centBounds.extend_atomic(centBounds);
      this->geomBounds.extend_atomic(geomBounds);
    }
    
    ObjectPartitionUnaligned::TaskBinParallel::TaskBinParallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, BezierRefList& prims, 
							       const LinearSpace3fa& space, const BBox3fa& geomBounds, const BBox3fa& centBounds) 
      : space(space), iter(prims), geomBounds(geomBounds), centBounds(centBounds), mapping(centBounds,space)
    {
      /* parallel binning */
      size_t numTasks = min(maxTasks,threadCount);
      scheduler->dispatchTask(threadIndex,numTasks,_task_bin_parallel,this,numTasks,"build::task_bin_parallel");
      
      /* reduction of bin information */
      BinInfo bins = binners[0];
      for (size_t i=1; i<numTasks; i++)
	bins.merge(binners[i]);
      
      /* calculation of best split */
      split = bins.best(prims,mapping);
    }
    
    void ObjectPartitionUnaligned::TaskBinParallel::task_bin_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      while (BezierRefList::item* block = iter.next())
	binners[taskIndex].bin(block->base(),block->size(),mapping);
    }
    
    template<>
    const ObjectPartitionUnaligned::Split ObjectPartitionUnaligned::find<true>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, BezierRefList& prims, const LinearSpace3fa& space, const PrimInfo& pinfo) 
    {
      /*! then perform parallel binning */
      return TaskBinParallel(threadIndex,threadCount,scheduler,prims,space,pinfo.geomBounds,pinfo.centBounds).split;
    }
    
    template<>
    void ObjectPartitionUnaligned::Split::split<false>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, 
						       PrimRefBlockAlloc<BezierPrim>& alloc, 
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
	  const Vec3fa center = prim.center(mapping.space);
	  const ssei bin = mapping.bin_unsafe(center);
	  
	  if (bin[dim] < pos) 
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
    
    ObjectPartitionUnaligned::TaskSplitParallel::TaskSplitParallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, const Split* split, PrimRefBlockAlloc<BezierPrim>& alloc, BezierRefList& prims, 
								   BezierRefList& lprims_o, PrimInfo& linfo_o, BezierRefList& rprims_o, PrimInfo& rinfo_o)
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
    
    void ObjectPartitionUnaligned::TaskSplitParallel::task_split_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      split->split<false>(threadIndex,threadCount,NULL,alloc,prims,lprims_o,linfos[taskIndex],rprims_o,rinfos[taskIndex]);
    }
    
    template<>
    void ObjectPartitionUnaligned::Split::split<true>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, 
						      PrimRefBlockAlloc<BezierPrim>& alloc, BezierRefList& prims, 
						      BezierRefList& lprims_o, PrimInfo& linfo_o, 
						      BezierRefList& rprims_o, PrimInfo& rinfo_o) const
    {
      TaskSplitParallel(threadIndex,threadCount,scheduler,this,alloc,prims,lprims_o,linfo_o,rprims_o,rinfo_o);
    }
  }
}
