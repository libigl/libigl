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

#include "heuristic_object_partition.h"

namespace embree
{
  namespace isa
  {
    //////////////////////////////////////////////////////////////////////////////
    //                        Bin Mapping                                       //
    //////////////////////////////////////////////////////////////////////////////
    
    __forceinline ObjectPartition::Mapping::Mapping(const PrimInfo& pinfo) 
    {
      num = min(maxBins,size_t(4.0f + 0.05f*pinfo.size()));
      const ssef diag = (ssef) pinfo.centBounds.size();
      scale = select(diag > ssef(1E-19f),rcp(diag) * ssef(0.99f*num),ssef(0.0f));
      ofs  = (ssef) pinfo.centBounds.lower;
    }
    
    __forceinline Vec3ia ObjectPartition::Mapping::bin(const Vec3fa& p) const 
    {
      const ssei i = floori((ssef(p)-ofs)*scale);
#if 1
      assert(i[0] >=0 && i[0] < num); 
      assert(i[1] >=0 && i[1] < num);
      assert(i[2] >=0 && i[2] < num);
      return Vec3ia(i);
#else
      return Vec3ia(clamp(i,ssei(0),ssei(num-1)));
#endif
    }

    __forceinline Vec3ia ObjectPartition::Mapping::bin_unsafe(const Vec3fa& p) const {
      return Vec3ia(floori((ssef(p)-ofs)*scale));
    }
    
    __forceinline bool ObjectPartition::Mapping::invalid(const int dim) const {
      return scale[dim] == 0.0f;
    }
    
    //////////////////////////////////////////////////////////////////////////////
    //                             Binning                                      //
    //////////////////////////////////////////////////////////////////////////////
    
    ObjectPartition::BinInfo::BinInfo() {
      clear();
    }

    __forceinline void ObjectPartition::BinInfo::clear() 
    {
      for (size_t i=0; i<maxBins; i++) {
	bounds[i][0] = bounds[i][1] = bounds[i][2] = bounds[i][3] = empty;
	counts[i] = 0;
      }
    }
    
    __forceinline void ObjectPartition::BinInfo::bin (const BezierPrim* prims, size_t N, const Mapping& mapping)
    {
      for (size_t i=0; i<N; i++)
      {
	const BBox3fa cbounds = prims[i].bounds();
	const Vec3fa  center  = prims[i].center();
	const ssei bin = ssei(mapping.bin(center));
	const int b0 = bin[0]; counts[b0][0]++; bounds[b0][0].extend(cbounds);
	const int b1 = bin[1]; counts[b1][1]++; bounds[b1][1].extend(cbounds);
	const int b2 = bin[2]; counts[b2][2]++; bounds[b2][2].extend(cbounds);
      }
    }
    
    __forceinline void ObjectPartition::BinInfo::bin (const PrimRef* prims, size_t num, const Mapping& mapping)
    {
      if (num == 0) return;
      
      size_t i; 
      for (i=0; i<num-1; i+=2)
      {
	/*! map even and odd primitive to bin */
	const BBox3fa prim0 = prims[i+0].bounds(); const Vec3fa center0 = Vec3fa(center2(prim0)); const Vec3ia bin0 = mapping.bin(center0); 
	const BBox3fa prim1 = prims[i+1].bounds(); const Vec3fa center1 = Vec3fa(center2(prim1)); const Vec3ia bin1 = mapping.bin(center1); 
	
	/*! increase bounds for bins for even primitive */
	const int b00 = bin0.x; counts[b00][0]++; bounds[b00][0].extend(prim0);
	const int b01 = bin0.y; counts[b01][1]++; bounds[b01][1].extend(prim0);
	const int b02 = bin0.z; counts[b02][2]++; bounds[b02][2].extend(prim0);
	
	/*! increase bounds of bins for odd primitive */
	const int b10 = bin1.x; counts[b10][0]++; bounds[b10][0].extend(prim1);
	const int b11 = bin1.y; counts[b11][1]++; bounds[b11][1].extend(prim1);
	const int b12 = bin1.z; counts[b12][2]++; bounds[b12][2].extend(prim1);
      }
      
      /*! for uneven number of primitives */
      if (i < num)
      {
	/*! map primitive to bin */
	const BBox3fa prim0 = prims[i].bounds(); const Vec3fa center0 = Vec3fa(center2(prim0)); const Vec3ia bin0 = mapping.bin(center0); 
	
	/*! increase bounds of bins */
	const int b00 = bin0.x; counts[b00][0]++; bounds[b00][0].extend(prim0);
	const int b01 = bin0.y; counts[b01][1]++; bounds[b01][1].extend(prim0);
	const int b02 = bin0.z; counts[b02][2]++; bounds[b02][2].extend(prim0);
      }
    }
    
    __forceinline void ObjectPartition::BinInfo::bin_copy (const PrimRef* prims, size_t num, const Mapping& mapping, PrimRef* dest)
    {
      if (num == 0) return;
      
      size_t i; 
      for (i=0; i<num-1; i+=2)
      {
	/*! map even and odd primitive to bin */
	const BBox3fa prim0 = prims[i+0].bounds(); const Vec3fa center0 = Vec3fa(center2(prim0)); const Vec3ia bin0 = mapping.bin(center0); 
	const BBox3fa prim1 = prims[i+1].bounds(); const Vec3fa center1 = Vec3fa(center2(prim1)); const Vec3ia bin1 = mapping.bin(center1); 
	
	/*! increase bounds for bins for even primitive */
	const int b00 = bin0.x; counts[b00][0]++; bounds[b00][0].extend(prim0);
	const int b01 = bin0.y; counts[b01][1]++; bounds[b01][1].extend(prim0);
	const int b02 = bin0.z; counts[b02][2]++; bounds[b02][2].extend(prim0);
	
	/*! increase bounds of bins for odd primitive */
	const int b10 = bin1.x; counts[b10][0]++; bounds[b10][0].extend(prim1);
	const int b11 = bin1.y; counts[b11][1]++; bounds[b11][1].extend(prim1);
	const int b12 = bin1.z; counts[b12][2]++; bounds[b12][2].extend(prim1);
	
	/*! copy to destination */
	dest[i+0] = prims[i+0];
	dest[i+1] = prims[i+1];
      }
      
      /*! for uneven number of primitives */
      if (i < num)
      {
	/*! map primitive to bin */
	const BBox3fa prim0 = prims[i].bounds(); const Vec3fa center0 = Vec3fa(center2(prim0)); const Vec3ia bin0 = mapping.bin(center0); 
	
	/*! increase bounds of bins */
	const int b00 = bin0.x; counts[b00][0]++; bounds[b00][0].extend(prim0);
	const int b01 = bin0.y; counts[b01][1]++; bounds[b01][1].extend(prim0);
	const int b02 = bin0.z; counts[b02][2]++; bounds[b02][2].extend(prim0);
	
	/*! copy to destination */
	dest[i+0] = prims[i+0];
      }
    }

    __forceinline void ObjectPartition::BinInfo::bin_copy (const PrimRef* prims, size_t begin, size_t end, const Mapping& mapping, PrimRef* dest)
    {
      bin_copy(prims+begin,end-begin,mapping,dest+begin);
    }

    __forceinline void ObjectPartition::BinInfo::bin(BezierRefList& prims, const Mapping& mapping)
    {
      BezierRefList::iterator i=prims;
      while (BezierRefList::item* block = i.next())
	bin(block->base(),block->size(),mapping);
    }
    
    __forceinline void ObjectPartition::BinInfo::bin(PrimRefList& prims, const Mapping& mapping)
    {
      PrimRefList::iterator i=prims;
      while (PrimRefList::item* block = i.next())
	bin(block->base(),block->size(),mapping);
    }
    
    __forceinline void ObjectPartition::BinInfo::merge (const BinInfo& other) 
    {
      for (size_t i=0; i<maxBins; i++) // FIXME: dont iterate over all bins
      {
	counts[i] += other.counts[i];
	bounds[i][0].extend(other.bounds[i][0]);
	bounds[i][1].extend(other.bounds[i][1]);
	bounds[i][2].extend(other.bounds[i][2]);
      }
    }

    __forceinline void ObjectPartition::BinInfo::merge (const BinInfo& other, size_t numBins)
    {
      for (size_t i=0; i<numBins; i++) 
      {
	counts[i] += other.counts[i];
	bounds[i][0].extend(other.bounds[i][0]);
	bounds[i][1].extend(other.bounds[i][1]);
	bounds[i][2].extend(other.bounds[i][2]);
      }
    }

    void ObjectPartition::BinInfo::reduce(const BinInfo binners[], size_t num, BinInfo& binner_o)
    {
      binner_o = binners[0];
      for (size_t tid=1; tid<num; tid++) 
      {
        const BinInfo& binner = binners[tid];
        for (size_t bin=0; bin<maxBins; bin++) 
        {
          binner_o.bounds[bin][0].extend(binner.bounds[bin][0]);
          binner_o.bounds[bin][1].extend(binner.bounds[bin][1]);
          binner_o.bounds[bin][2].extend(binner.bounds[bin][2]);
          binner_o.counts[bin] += binner.counts[bin];
        }
      }
    }

    void ObjectPartition::BinInfo::reduce2(const BinInfo binners[], size_t num, BinInfo& binner_o)
    {
      binner_o = binners[0];
      for (size_t tid=1; tid<num; tid++) 
      {
        const BinInfo& binner = binners[tid];
        for (size_t bin=0; bin<16; bin++) 
        {
          binner_o.bounds[bin][0].extend(binner.bounds[bin][0]);
          binner_o.bounds[bin][1].extend(binner.bounds[bin][1]);
          binner_o.bounds[bin][2].extend(binner.bounds[bin][2]);
          binner_o.counts[bin] += binner.counts[bin];
        }
      }
    }
    
    __forceinline ObjectPartition::Split ObjectPartition::BinInfo::best(const Mapping& mapping, const size_t blocks_shift)
    {
      /* sweep from right to left and compute parallel prefix of merged bounds */
      ssef rAreas[maxBins];
      ssei rCounts[maxBins];
      ssei count = 0; BBox3fa bx = empty; BBox3fa by = empty; BBox3fa bz = empty;
      for (size_t i=mapping.size()-1; i>0; i--)
      {
	count += counts[i];
	rCounts[i] = count;
	bx.extend(bounds[i][0]); rAreas[i][0] = halfArea(bx);
	by.extend(bounds[i][1]); rAreas[i][1] = halfArea(by);
	bz.extend(bounds[i][2]); rAreas[i][2] = halfArea(bz);
      }
      
      /* sweep from left to right and compute SAH */
      ssei blocks_add = (1 << blocks_shift)-1;
      ssei ii = 1; ssef vbestSAH = pos_inf; ssei vbestPos = 0; 
      count = 0; bx = empty; by = empty; bz = empty;
      for (size_t i=1; i<mapping.size(); i++, ii+=1)
      {
	count += counts[i-1];
	bx.extend(bounds[i-1][0]); float Ax = halfArea(bx);
	by.extend(bounds[i-1][1]); float Ay = halfArea(by);
	bz.extend(bounds[i-1][2]); float Az = halfArea(bz);
	const ssef lArea = ssef(Ax,Ay,Az,Az);
	const ssef rArea = rAreas[i];
	const ssei lCount = (count     +blocks_add) >> blocks_shift;
	const ssei rCount = (rCounts[i]+blocks_add) >> blocks_shift;
	const ssef sah = lArea*ssef(lCount) + rArea*ssef(rCount);
	vbestPos = select(sah < vbestSAH,ii ,vbestPos);
	vbestSAH = select(sah < vbestSAH,sah,vbestSAH);
      }
      
      /* find best dimension */
      float bestSAH = inf;
      int   bestDim = -1;
      int   bestPos = 0;
      int   bestLeft = 0;
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
      
      return ObjectPartition::Split(bestSAH,bestDim,bestPos,mapping);
    }

    template<>
    const ObjectPartition::Split ObjectPartition::find<false>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, BezierRefList& prims, const PrimInfo& pinfo, const size_t logBlockSize)
    {
      BinInfo binner;
      const Mapping mapping(pinfo);
      binner.bin(prims,mapping);
      return binner.best(mapping,logBlockSize);
    }
    
    template<>
    const ObjectPartition::Split ObjectPartition::find<false>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefList& prims, const PrimInfo& pinfo, const size_t logBlockSize)
    {
      BinInfo binner;
      const Mapping mapping(pinfo);
      binner.bin(prims,mapping);
      return binner.best(mapping,logBlockSize);
    }
    
    const ObjectPartition::Split ObjectPartition::find(PrimRef *__restrict__ const prims, const size_t begin, const size_t end, const PrimInfo& pinfo, const size_t logBlockSize)
    {
      BinInfo binner;
      const Mapping mapping(pinfo);
      binner.bin(prims+begin,end-begin,mapping);
      return binner.best(mapping,logBlockSize);
    }

    template<>
    const ObjectPartition::Split ObjectPartition::find<false>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefList& prims, const PrimInfo& pinfo, const size_t logBlockSize, SplitInfo& sinfo_o)
    {
      BinInfo binner;
      const Mapping mapping(pinfo);
      binner.bin(prims,mapping);
      const ObjectPartition::Split split = binner.best(mapping,logBlockSize);
      binner.getSplitInfo(mapping,split,sinfo_o);
      return split;
    }

    template<>
    const std::pair<BBox3fa,BBox3fa> ObjectPartition::computePrimInfoMB<false>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, Scene* scene, BezierRefList& prims)
    {
      BBox3fa bounds0 = empty;
      BBox3fa bounds1 = empty;
      for (BezierRefList::block_iterator_unsafe i = prims; i; i++) 
      {
        const BezierCurves* curves = scene->getBezierCurves(i->geomID<0>());
        bounds0.extend(curves->bounds(i->primID<0>(),0));
        bounds1.extend(curves->bounds(i->primID<0>(),1));
      }
      return std::pair<BBox3fa,BBox3fa>(bounds0,bounds1);
    }
    
    template<>
    const std::pair<BBox3fa,BBox3fa> ObjectPartition::computePrimInfoMB<true>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, Scene* scene, BezierRefList& prims)
    {
      const TaskPrimInfoMBParallel bounds(threadIndex,threadCount,scheduler,scene,prims);
      return std::pair<BBox3fa,BBox3fa>(bounds.bounds0,bounds.bounds1);
    }

    ObjectPartition::TaskPrimInfoMBParallel::TaskPrimInfoMBParallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, Scene* scene, BezierRefList& prims) 
      : scene(scene), iter(prims), bounds0(empty), bounds1(empty)
    {
      size_t numTasks = min(maxTasks,threadCount);
      scheduler->dispatchTask(threadIndex,numTasks,_task_bound_parallel,this,numTasks,"build::task_bound_parallel");
    }
    
    void ObjectPartition::TaskPrimInfoMBParallel::task_bound_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      size_t N = 0;
      BBox3fa bounds0 = empty;
      BBox3fa bounds1 = empty;
      while (BezierRefList::item* block = iter.next()) 
      {
	for (size_t i=0; i<block->size(); i++) 
        {
          const BezierPrim& ref = block->at(i);
          const BezierCurves* curves = scene->getBezierCurves(ref.geomID<0>());
          bounds0.extend(curves->bounds(ref.primID<0>(),0));
          bounds1.extend(curves->bounds(ref.primID<0>(),1));
	}
      }
      this->bounds0.extend_atomic(bounds0);
      this->bounds1.extend_atomic(bounds1);
    }
    
    //////////////////////////////////////////////////////////////////////////////
    //                         Parallel Binning                                 //
    //////////////////////////////////////////////////////////////////////////////
    
    template<typename List>
    ObjectPartition::TaskBinParallel<List>::TaskBinParallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, List& prims, const PrimInfo& pinfo, const size_t logBlockSize) 
      : iter(prims)
    {
      /* parallel binning */			
      size_t numTasks = min(maxTasks,threadCount);
      new (&mapping) Mapping(pinfo);
      scheduler->dispatchTask(threadIndex,numTasks,_task_bin_parallel,this,numTasks,"build::task_bin_parallel");
      
      /* reduction of bin informations */
      binner = binners[0];
      for (size_t i=1; i<numTasks; i++)
	binner.merge(binners[i]);
      
      /* calculation of best split */
      split = binner.best(mapping,logBlockSize);
    }
    
    template<typename List>
    void ObjectPartition::TaskBinParallel<List>::task_bin_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      while (typename List::item* block = iter.next())
	binners[taskIndex].bin(block->base(),block->size(),mapping);
    }
    
    template<>
    const ObjectPartition::Split ObjectPartition::find<true>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, BezierRefList& prims, const PrimInfo& pinfo, const size_t logBlockSize) {
      return TaskBinParallel<BezierRefList>(threadIndex,threadCount,scheduler,prims,pinfo,logBlockSize).split;
    }
    
    template<>
    const ObjectPartition::Split ObjectPartition::find<true>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefList& prims, const PrimInfo& pinfo, const size_t logBlockSize) {
      return TaskBinParallel<PrimRefList>(threadIndex,threadCount,scheduler,prims,pinfo,logBlockSize).split;
    }
    
    template<>
    const ObjectPartition::Split ObjectPartition::find<true>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefList& prims, const PrimInfo& pinfo, const size_t logBlockSize, SplitInfo& sinfo_o) 
    {
      TaskBinParallel<PrimRefList> task(threadIndex,threadCount,scheduler,prims,pinfo,logBlockSize);
      task.binner.getSplitInfo(task.mapping,task.split,sinfo_o);
      return task.split;
    }
    

    //////////////////////////////////////////////////////////////////////////////
    //                             Splitting                                    //
    //////////////////////////////////////////////////////////////////////////////
    
    template<>
    void ObjectPartition::Split::split<false>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, 
					      PrimRefBlockAlloc<BezierPrim>& alloc, 
					      BezierRefList& prims, 
					      BezierRefList& lprims_o, PrimInfo& linfo_o, 
					      BezierRefList& rprims_o, PrimInfo& rinfo_o) const
    {
      assert(valid());
      BezierRefList::item* lblock = lprims_o.insert(alloc.malloc(threadIndex));
      BezierRefList::item* rblock = rprims_o.insert(alloc.malloc(threadIndex));
      linfo_o.reset();
      rinfo_o.reset();

      while (BezierRefList::item* block = prims.take()) 
      {
	for (size_t i=0; i<block->size(); i++) 
	{
	  const BezierPrim& prim = block->at(i); 
	  const Vec3fa center = prim.center2();
	  const ssei bin = ssei(mapping.bin_unsafe(center));
	  
	  if (bin[dim] < pos) 
	  {
	    linfo_o.add(prim.bounds(),center);
	    if (likely(lblock->insert(prim))) continue; 
	    lblock = lprims_o.insert(alloc.malloc(threadIndex));
	    lblock->insert(prim);
	  } 
	  else 
	  {
	    rinfo_o.add(prim.bounds(),center);
	    if (likely(rblock->insert(prim))) continue;
	    rblock = rprims_o.insert(alloc.malloc(threadIndex));
	    rblock->insert(prim);
	  }
	}
	alloc.free(threadIndex,block);
      }
    }
    
    template<>
    void ObjectPartition::Split::split<false>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, 
					      PrimRefBlockAlloc<PrimRef>& alloc, 
					      PrimRefList& prims, 
					      PrimRefList& lprims_o, PrimInfo& linfo_o, 
					      PrimRefList& rprims_o, PrimInfo& rinfo_o) const
    {
      assert(valid());
      PrimRefList::item* lblock = lprims_o.insert(alloc.malloc(threadIndex));
      PrimRefList::item* rblock = rprims_o.insert(alloc.malloc(threadIndex));
      linfo_o.reset();
      rinfo_o.reset();

      size_t numLeft = 0; CentGeomBBox3fa leftBounds(empty);
      size_t numRight = 0; CentGeomBBox3fa rightBounds(empty);
      
      while (PrimRefList::item* block = prims.take()) 
      {
	for (size_t i=0; i<block->size(); i++) 
	{
	  const PrimRef& prim = block->at(i); 
	  const Vec3fa center = center2(prim.bounds());
	  const ssei bin = ssei(mapping.bin_unsafe(center));

	  if (bin[dim] < pos) 
	  {
	    leftBounds.extend(prim.bounds()); numLeft++;
	    //linfo_o.add(prim.bounds(),center);
	    //if (++lblock->num > PrimRefBlock::blockSize)
	    //lblock = lprims_o.insert(alloc.malloc(threadIndex));
	    if (likely(lblock->insert(prim))) continue; 
	    lblock = lprims_o.insert(alloc.malloc(threadIndex));
	    lblock->insert(prim);
	  } 
	  else 
	  {
	    rightBounds.extend(prim.bounds()); numRight++;
	    //rinfo_o.add(prim.bounds(),center);
	    //if (++rblock->num > PrimRefBlock::blockSize)
	    //rblock = rprims_o.insert(alloc.malloc(threadIndex));
	    if (likely(rblock->insert(prim))) continue;
	    rblock = rprims_o.insert(alloc.malloc(threadIndex));
	    rblock->insert(prim);
	  }
	}
	alloc.free(threadIndex,block);
      }

      linfo_o.add(leftBounds.geomBounds,leftBounds.centBounds,numLeft);
      rinfo_o.add(rightBounds.geomBounds,rightBounds.centBounds,numRight);
    }
        
    void ObjectPartition::Split::partition(PrimRef *__restrict__ const prims, const size_t begin, const size_t end, PrimInfo& left, PrimInfo& right) const
    {
      assert(valid());
      CentGeomBBox3fa local_left(empty);
      CentGeomBBox3fa local_right(empty);
      
      assert(begin <= end);
      PrimRef* l = prims + begin;
      PrimRef* r = prims + end - 1;

      while(1)
      {
	while (likely(l <= r && mapping.bin_unsafe(center2(l->bounds()))[dim] < pos)) 
	{
	  local_left.extend(l->bounds());
	  ++l;
	}
	while (likely(l <= r && mapping.bin_unsafe(center2(r->bounds()))[dim] >= pos)) 
	{
	  local_right.extend(r->bounds());
	  --r;
	}
	if (r<l) break;
	
	const BBox3fa bl = l->bounds();
	const BBox3fa br = r->bounds();
	local_left.extend(br);
	local_right.extend(bl);
	*(BBox3fa*)l = br;
	*(BBox3fa*)r = bl;
	l++; r--;
      }
      
      unsigned int center = l - prims;
      new (&left ) PrimInfo(begin,center,local_left.geomBounds,local_left.centBounds);
      new (&right) PrimInfo(center,end,local_right.geomBounds,local_right.centBounds);
      assert(area(left.geomBounds) >= 0.0f);
      assert(area(right.geomBounds) >= 0.0f);
    }
    
    template<typename Prim>
    ObjectPartition::TaskSplitParallel<Prim>::TaskSplitParallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, const Split* split, PrimRefBlockAlloc<Prim>& alloc, List& prims, 
								List& lprims_o, PrimInfo& linfo_o, List& rprims_o, PrimInfo& rinfo_o)
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
    
    template<typename Prim>
    void ObjectPartition::TaskSplitParallel<Prim>::task_split_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      split->split<false>(threadIndex,threadCount,NULL,alloc,prims,lprims_o,linfos[taskIndex],rprims_o,rinfos[taskIndex]);
    }
    
    template<>
    void ObjectPartition::Split::split<true>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, 
					     PrimRefBlockAlloc<BezierPrim>& alloc, BezierRefList& prims, 
					     BezierRefList& lprims_o, PrimInfo& linfo_o, 
					     BezierRefList& rprims_o, PrimInfo& rinfo_o) const
    {
      TaskSplitParallel<BezierPrim>(threadIndex,threadCount,scheduler,this,alloc,prims,lprims_o,linfo_o,rprims_o,rinfo_o);
    }
    
    template<>
    void ObjectPartition::Split::split<true>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, 
					     PrimRefBlockAlloc<PrimRef>& alloc, PrimRefList& prims, 
					     PrimRefList& lprims_o, PrimInfo& linfo_o, 
					     PrimRefList& rprims_o, PrimInfo& rinfo_o) const
    {
      TaskSplitParallel<PrimRef>(threadIndex,threadCount,scheduler,this,alloc,prims,lprims_o,linfo_o,rprims_o,rinfo_o);
    }
    
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    
    void ObjectPartition::ParallelBinner::parallelBinning(size_t threadID, size_t numThreads)
    {
      BinInfo& bin16 = global_bin16[threadID];
      const size_t startID = pinfo.begin + (threadID+0)*pinfo.size()/numThreads;
      const size_t endID   = pinfo.begin + (threadID+1)*pinfo.size()/numThreads;
      bin16.clear();
      bin16.bin_copy(src,startID,endID,mapping,dst);
    }
    
    float ObjectPartition::ParallelBinner::find(const PrimInfo& pinfo, const PrimRef* src, PrimRef* dst, const size_t logBlockSize, const size_t threadID, const size_t numThreads, LockStepTaskScheduler* scheduler) 
    {
      this->pinfo = pinfo;
      mapping = Mapping(pinfo);
      left.reset();
      right.reset();
      this->src = src;
      this->dst = dst;
      scheduler->dispatchTask(task_parallelBinning, this, threadID, numThreads );
      
      /* reduce binning information from all threads */
      bin16 = global_bin16[0];
      for (size_t i=1; i<numThreads; i++)
	bin16.merge(global_bin16[i]);

      split = bin16.best(mapping,logBlockSize);
      return split.sah;
    }
    
    void ObjectPartition::ParallelBinner::parallelPartition(size_t threadID, size_t numThreads)
    {
      const size_t startID = pinfo.begin + (threadID+0)*pinfo.size()/numThreads;
      const size_t endID   = pinfo.begin + (threadID+1)*pinfo.size()/numThreads;
      
      /* load binning function */
      const unsigned int splitPos = split.pos;
      const unsigned int splitDim = split.dim;
      const float centroidBase = mapping.ofs[splitDim];
      const float centroidScale = mapping.scale[splitDim];
      
      /* compute items per thread that go to the 'left' and to the 'right' */
      int lnum[maxBins];
      lnum[0] = global_bin16[threadID].counts[0][splitDim];
      for (size_t i=1; i<mapping.size(); i++)
        lnum[i] = lnum[i-1] + global_bin16[threadID].counts[i][splitDim];

      size_t numLeft = bin16.getNumLeft(split);
      
      const unsigned int localNumLeft = lnum[splitPos-1];
      const unsigned int localNumRight = (endID-startID) - localNumLeft;
      
      const unsigned int startLeft  = lCounter.add(localNumLeft);
      const unsigned int startRight = rCounter.add(localNumRight);
      
      PrimRef* __restrict__ src = (PrimRef*)this->src;
      PrimRef* __restrict__ dstLeft = dst + pinfo.begin + startLeft;
      PrimRef* __restrict__ dstRight = dst + pinfo.begin + startRight + numLeft;
      
      /* split into left and right */
      CentGeomBBox3fa leftBounds(empty);
      CentGeomBBox3fa rightBounds(empty);
      
      for (size_t i=startID; i<endID; i++)
      {
        if (likely(mapping.bin_unsafe(center2(src[i].bounds()))[splitDim] < splitPos)) {
          leftBounds.extend(src[i].bounds()); 
          *dstLeft++ = src[i];
        } else {
          rightBounds.extend(src[i].bounds()); 
          *dstRight++ = src[i];
        }
      }
      
      left .extend_atomic(leftBounds); 
      right.extend_atomic(rightBounds);  
    }
    
    void ObjectPartition::ParallelBinner::partition(const PrimInfo& pinfo, const PrimRef* src, PrimRef* dst, PrimInfo& leftChild, PrimInfo& rightChild, const size_t threadID, const size_t numThreads, LockStepTaskScheduler* scheduler)
    {
      left.reset(); lCounter.reset(0);
      right.reset(); rCounter.reset(0); 
      this->src = src;
      this->dst = dst;
      scheduler->dispatchTask(task_parallelPartition, this, threadID, numThreads);
      size_t numLeft = bin16.getNumLeft(split);
      unsigned center = pinfo.begin + numLeft;
      assert(lCounter == numLeft);
      assert(rCounter == pinfo.size() - lCounter);
      new (&leftChild ) PrimInfo(pinfo.begin,center,left.geomBounds,left.centBounds);
      new (&rightChild) PrimInfo(center,pinfo.end,right.geomBounds,right.centBounds);
    }
  }
}
