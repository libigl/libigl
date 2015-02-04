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

#include "heuristic_spatial_split.h"

namespace embree
{
  namespace isa
  {
    __forceinline void splitTriangle(const PrimRef& prim, int dim, float pos, 
				     const Vec3fa& a, const Vec3fa& b, const Vec3fa& c, PrimRef& left_o, PrimRef& right_o)
    {
      BBox3fa left = empty, right = empty;
      const Vec3fa v[3] = { a,b,c };
      
      /* clip triangle to left and right box by processing all edges */
      Vec3fa v1 = v[2];
      for (size_t i=0; i<3; i++)
      {
	Vec3fa v0 = v1; v1 = v[i];
	float v0d = v0[dim], v1d = v1[dim];
	
	if (v0d <= pos) left. extend(v0); // this point is on left side
	if (v0d >= pos) right.extend(v0); // this point is on right side
	
	if ((v0d < pos && pos < v1d) || (v1d < pos && pos < v0d)) // the edge crosses the splitting location
	{
	  assert((v1d-v0d) != 0.0f);
	  Vec3fa c = v0 + (pos-v0d)/(v1d-v0d)*(v1-v0);
	  left.extend(c);
	  right.extend(c);
	}
      }
      //assert(!left.empty());  // happens if split does not hit triangle
      //assert(!right.empty()); // happens if split does not hit triangle
      
      /* clip against current bounds */
      BBox3fa bounds = prim.bounds();
      BBox3fa cleft (max(left .lower,bounds.lower),min(left .upper,bounds.upper));
      BBox3fa cright(max(right.lower,bounds.lower),min(right.upper,bounds.upper));
      
      new (&left_o ) PrimRef(cleft, prim.geomID(), prim.primID());
      new (&right_o) PrimRef(cright,prim.geomID(), prim.primID());
    }
    
    __forceinline SpatialSplit::Mapping::Mapping(const PrimInfo& pinfo) 
    {
      const ssef lower = (ssef) pinfo.geomBounds.lower;
      const ssef upper = (ssef) pinfo.geomBounds.upper;
      const sseb ulpsized = upper - lower <= max(ssef(1E-19f),128.0f*ssef(ulp)*max(abs(lower),abs(upper)));
      const ssef diag = (ssef) pinfo.geomBounds.size();
      scale = select(ulpsized,ssef(0.0f),rcp(diag) * ssef(BINS * 0.99f));
      ofs  = (ssef) pinfo.geomBounds.lower;
    }
    
    __forceinline ssei SpatialSplit::Mapping::bin(const Vec3fa& p) const 
    {
      const ssei i = floori((ssef(p)-ofs)*scale);
      return clamp(i,ssei(0),ssei(BINS-1));
    }
    
    __forceinline float SpatialSplit::Mapping::pos(const int bin, const int dim) const {
      return float(bin)/scale[dim]+ofs[dim];
    }
    
    __forceinline bool SpatialSplit::Mapping::invalid(const int dim) const {
      return scale[dim] == 0.0f;
    }
    
    __forceinline SpatialSplit::BinInfo::BinInfo()
    {
      for (size_t i=0; i<BINS; i++) {
	bounds[i][0] = bounds[i][1] = bounds[i][2] = bounds[i][3] = empty;
	numBegin[i] = numEnd[i] = 0;
      }
    }
    
    __forceinline void SpatialSplit::BinInfo::bin (Scene* scene, const BezierPrim* prims, size_t N, const PrimInfo& pinfo, const Mapping& mapping)
    {
      for (size_t i=0; i<N; i++)
      {
	const ssei bin0 = mapping.bin(min(prims[i].p0,prims[i].p3));
	const ssei bin1 = mapping.bin(max(prims[i].p0,prims[i].p3));
	
	for (size_t dim=0; dim<3; dim++) 
	{
	  size_t bin;
	  BezierPrim curve = prims[i];
	  for (bin=bin0[dim]; bin<bin1[dim]; bin++)
	  {
	    const float pos = mapping.pos(bin+1,dim);
	    BezierPrim bincurve,restcurve; 
	    if (curve.split(dim,pos,bincurve,restcurve)) {
	      bounds[bin][dim].extend(bincurve.bounds());
	      curve = restcurve;
	    }
	  }
	  numBegin[bin0[dim]][dim]++;
	  numEnd  [bin1[dim]][dim]++;
	  bounds  [bin][dim].extend(curve.bounds());
	}
      }
    }

    __forceinline void SpatialSplit::BinInfo::bin(Scene* scene, const PrimRef* prims, size_t N, const PrimInfo& pinfo, const Mapping& mapping)
    {
      for (size_t i=0; i<N; i++)
      {
	const PrimRef prim = prims[i];
	TriangleMesh* mesh = (TriangleMesh*) scene->get(prim.geomID());
	TriangleMesh::Triangle tri = mesh->triangle(prim.primID());
	const Vec3fa v0 = mesh->vertex(tri.v[0]);
	const Vec3fa v1 = mesh->vertex(tri.v[1]);
	const Vec3fa v2 = mesh->vertex(tri.v[2]);
	const ssei bin0 = mapping.bin(prim.bounds().lower);
	const ssei bin1 = mapping.bin(prim.bounds().upper);

	for (size_t dim=0; dim<3; dim++) 
	{
	  size_t bin;
	  PrimRef rest = prim;
	  size_t l = bin0[dim];
	  size_t r = bin1[dim];
	  for (bin=bin0[dim]; bin<bin1[dim]; bin++) 
	  {
	    const float pos = mapping.pos(bin+1,dim);
	    
	    PrimRef left,right;
	    splitTriangle(rest,dim,pos,v0,v1,v2,left,right);
	    if (left.bounds().empty()) l++;
	    
	    bounds[bin][dim].extend(left.bounds());
	    rest = right;
	  }
	  if (rest.bounds().empty()) r--;
	  //numBegin[bin0[dim]][dim]++;
	  //numEnd  [bin1[dim]][dim]++;
	  numBegin[l][dim]++;
	  numEnd  [r][dim]++;
	  bounds  [bin][dim].extend(rest.bounds());
	}
      }
    }
    
    __forceinline void SpatialSplit::BinInfo::bin(Scene* scene, BezierRefList& prims, const PrimInfo& pinfo, const Mapping& mapping)
    {
      BezierRefList::iterator i=prims;
      while (BezierRefList::item* block = i.next())
	bin(scene,block->base(),block->size(),pinfo,mapping);
    }

    __forceinline void SpatialSplit::BinInfo::bin(Scene* scene, TriRefList& prims, const PrimInfo& pinfo, const Mapping& mapping)
    {
      TriRefList::iterator i=prims;
      while (TriRefList::item* block = i.next())
	bin(scene,block->base(),block->size(),pinfo,mapping);
    }
    
    __forceinline void SpatialSplit::BinInfo::merge (const BinInfo& other) // FIXME: dont iterate over all bins
    {
      for (size_t i=0; i<BINS; i++) 
      {
	numBegin[i] += other.numBegin[i];
	numEnd  [i] += other.numEnd  [i];
	bounds[i][0].extend(other.bounds[i][0]);
	bounds[i][1].extend(other.bounds[i][1]);
	bounds[i][2].extend(other.bounds[i][2]);
      }
    }
    
    __forceinline SpatialSplit::Split SpatialSplit::BinInfo::best(const PrimInfo& pinfo, const Mapping& mapping, const size_t blocks_shift)
    {
      /* sweep from right to left and compute parallel prefix of merged bounds */
      ssef rAreas[BINS];
      ssei rCounts[BINS];
      ssei count = 0; BBox3fa bx = empty; BBox3fa by = empty; BBox3fa bz = empty;
      for (size_t i=BINS-1; i>0; i--)
      {
	count += numEnd[i];
	rCounts[i] = count;
	bx.extend(bounds[i][0]); rAreas[i][0] = halfArea(bx);
	by.extend(bounds[i][1]); rAreas[i][1] = halfArea(by);
	bz.extend(bounds[i][2]); rAreas[i][2] = halfArea(bz);
      }
      
      /* sweep from left to right and compute SAH */
      ssei blocks_add = (1 << blocks_shift)-1;
      ssei ii = 1; ssef vbestSAH = pos_inf; ssei vbestPos = 0;
      count = 0; bx = empty; by = empty; bz = empty;
      for (size_t i=1; i<BINS; i++, ii+=1)
      {
	count += numBegin[i-1];
	bx.extend(bounds[i-1][0]); float Ax = halfArea(bx);
	by.extend(bounds[i-1][1]); float Ay = halfArea(by);
	bz.extend(bounds[i-1][2]); float Az = halfArea(bz);
	const ssef lArea = ssef(Ax,Ay,Az,Az);
	const ssef rArea = rAreas[i];
	const ssei lCount = (count     +blocks_add) >> blocks_shift;
	const ssei rCount = (rCounts[i]+blocks_add) >> blocks_shift;
	const ssef sah = lArea*ssef(lCount) + rArea*ssef(rCount);
	vbestPos  = select(sah < vbestSAH,ii ,vbestPos);
	vbestSAH  = select(sah < vbestSAH,sah,vbestSAH);
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
      
      /* return invalid split if no split found */
      if (bestDim == -1) 
	return Split(inf,-1,0,mapping);
      
      /* return best found split */
      return Split(bestSAH,bestDim,bestPos,mapping);
    }
    
    template<>
    const SpatialSplit::Split SpatialSplit::find<false>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, Scene* scene, BezierRefList& prims, const PrimInfo& pinfo, const size_t logBlockSize)
    {
      BinInfo binner;
      Mapping mapping(pinfo);
      binner.bin(scene,prims,pinfo,mapping);
      return binner.best(pinfo,mapping,logBlockSize);
    }

    template<>
    const SpatialSplit::Split SpatialSplit::find<false>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, Scene* scene, TriRefList& prims, const PrimInfo& pinfo, const size_t logBlockSize)
    {
      BinInfo binner;
      Mapping mapping(pinfo);
      binner.bin(scene,prims,pinfo,mapping);
      return binner.best(pinfo,mapping,logBlockSize);
    }
    
    //////////////////////////////////////////////////////////////////////////////
    //                         Parallel Binning                                 //
    //////////////////////////////////////////////////////////////////////////////
   
    template<typename List>
    SpatialSplit::TaskBinParallel<List>::TaskBinParallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, Scene* scene, List& prims, const PrimInfo& pinfo, const Mapping& mapping, const size_t logBlockSize) 
      : scene(scene), iter(prims), pinfo(pinfo), mapping(mapping)
    {
      /* parallel binning */
      size_t numTasks = min(maxTasks,threadCount);
      scheduler->dispatchTask(threadIndex,numTasks,_task_bin_parallel,this,numTasks,"build::task_bin_parallel");
      
      /* reduction of bin informations */
      BinInfo bins = binners[0];
      for (size_t i=1; i<numTasks; i++)
	bins.merge(binners[i]);
      
      /* calculation of best split */
      split = bins.best(pinfo,mapping,logBlockSize);
    }
    
    template<typename List>
    void SpatialSplit::TaskBinParallel<List>::task_bin_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      while (typename List::item* block = iter.next())
	binners[taskIndex].bin(scene,block->base(),block->size(),pinfo,mapping);
    }
    
    template<>
    const SpatialSplit::Split SpatialSplit::find<true>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, Scene* scene, BezierRefList& prims, const PrimInfo& pinfo, const size_t logBlockSize) 
    {
      const Mapping mapping(pinfo);
      return TaskBinParallel<BezierRefList>(threadIndex,threadCount,scheduler,scene,prims,pinfo,mapping,logBlockSize).split;
    }

    template<>
    const SpatialSplit::Split SpatialSplit::find<true>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, Scene* scene, TriRefList& prims, const PrimInfo& pinfo, const size_t logBlockSize) 
    {
      const Mapping mapping(pinfo);
      return TaskBinParallel<TriRefList>(threadIndex,threadCount,scheduler,scene,prims,pinfo,mapping,logBlockSize).split;
    }
    
    //////////////////////////////////////////////////////////////////////////////
    //                             Splitting                                    //
    //////////////////////////////////////////////////////////////////////////////
    
    template<>
    void SpatialSplit::Split::split<false>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefBlockAlloc<BezierPrim>& alloc, 
					   Scene* scene, BezierRefList& prims, 
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
	  const int bin0 = mapping.bin(min(prim.p0,prim.p3))[dim];
	  const int bin1 = mapping.bin(max(prim.p0,prim.p3))[dim];
	  
	  /* sort to the left side */
	  if (bin0 < pos && bin1 < pos) // FIXME: optimize
	  {
	    linfo_o.add(prim.bounds(),prim.center());
	    if (likely(lblock->insert(prim))) continue; 
	    lblock = lprims_o.insert(alloc.malloc(threadIndex));
	    lblock->insert(prim);
	    continue;
	  }
	  
	  /* sort to the right side */
	  if (bin0 >= pos && bin1 >= pos)// FIXME: optimize
	  {
	    rinfo_o.add(prim.bounds(),prim.center());
	    if (likely(rblock->insert(prim))) continue;
	    rblock = rprims_o.insert(alloc.malloc(threadIndex));
	    rblock->insert(prim);
	    continue;
	  }
	  
	  /* split and sort to left and right */
	  BezierPrim left,right;
	  float fpos = mapping.pos(pos,dim);
	  if (prim.split(dim,fpos,left,right)) 
	  {
	    if (!left.bounds().empty()) {
	      linfo_o.add(left.bounds(),left.center());
	      if (!lblock->insert(left)) {
		lblock = lprims_o.insert(alloc.malloc(threadIndex));
		lblock->insert(left);
	      }
	    }
	    if (!right.bounds().empty()) {
	      rinfo_o.add(right.bounds(),right.center());
	      if (!rblock->insert(right)) {
		rblock = rprims_o.insert(alloc.malloc(threadIndex));
		rblock->insert(right);
	      }
	    }
	    continue;
	  }
	  
	  /* insert to left side as fallback */
	  linfo_o.add(prim.bounds(),prim.center());
	  if (!lblock->insert(prim)) {
	    lblock = lprims_o.insert(alloc.malloc(threadIndex));
	    lblock->insert(prim);
	  }
	}
	alloc.free(threadIndex,block);
      }
    }

    template<>
    void SpatialSplit::Split::split<false>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefBlockAlloc<PrimRef>& alloc, 
					   Scene* scene, TriRefList& prims, 
					   TriRefList& lprims_o, PrimInfo& linfo_o, 
					   TriRefList& rprims_o, PrimInfo& rinfo_o) const
    {
      assert(valid());
      TriRefList::item* lblock = lprims_o.insert(alloc.malloc(threadIndex));
      TriRefList::item* rblock = rprims_o.insert(alloc.malloc(threadIndex));
      linfo_o.reset();
      rinfo_o.reset();
    
      /* sort each primitive to left, right, or left and right */
      while (atomic_set<PrimRefBlock>::item* block = prims.take()) 
      {
	for (size_t i=0; i<block->size(); i++) 
	{
	  const PrimRef& prim = block->at(i); 
	  const BBox3fa bounds = prim.bounds();
	  const int bin0 = mapping.bin(bounds.lower)[dim];
	  const int bin1 = mapping.bin(bounds.upper)[dim];

	  /* sort to the left side */
	  if (bin1 < pos)
	  {
	    linfo_o.add(bounds,center2(bounds));
	    if (likely(lblock->insert(prim))) continue; 
	    lblock = lprims_o.insert(alloc.malloc(threadIndex));
	    lblock->insert(prim);
	    continue;
	  }
	  
	  /* sort to the right side */
	  if (bin0 >= pos)
	  {
	    rinfo_o.add(bounds,center2(bounds));
	    if (likely(rblock->insert(prim))) continue;
	    rblock = rprims_o.insert(alloc.malloc(threadIndex));
	    rblock->insert(prim);
	    continue;
	  }
	  
	  /* split and sort to left and right */
	  TriangleMesh* mesh = (TriangleMesh*) scene->get(prim.geomID());
	  TriangleMesh::Triangle tri = mesh->triangle(prim.primID());
	  const Vec3fa v0 = mesh->vertex(tri.v[0]);
	  const Vec3fa v1 = mesh->vertex(tri.v[1]);
	  const Vec3fa v2 = mesh->vertex(tri.v[2]);
	  
	  PrimRef left,right;
	  float fpos = mapping.pos(pos,dim);
	  splitTriangle(prim,dim,fpos,v0,v1,v2,left,right);
	
	  if (!left.bounds().empty()) {
	    linfo_o.add(left.bounds(),center2(left.bounds()));
	    if (!lblock->insert(left)) {
	      lblock = lprims_o.insert(alloc.malloc(threadIndex));
	      lblock->insert(left);
	    }
	  }
	  
	  if (!right.bounds().empty()) {
	    rinfo_o.add(right.bounds(),center2(right.bounds()));
	    if (!rblock->insert(right)) {
	      rblock = rprims_o.insert(alloc.malloc(threadIndex));
	      rblock->insert(right);
	    }
	  }
	}
	alloc.free(threadIndex,block);
      }
    }
    
    template<typename Prim>
    SpatialSplit::TaskSplitParallel<Prim>::TaskSplitParallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, const Split* split, 
							     PrimRefBlockAlloc<Prim>& alloc, Scene* scene, List& prims, 
							     List& lprims_o, PrimInfo& linfo_o, 
							     List& rprims_o, PrimInfo& rinfo_o)
      : split(split), alloc(alloc), scene(scene), prims(prims), lprims_o(lprims_o), linfo_o(linfo_o), rprims_o(rprims_o), rinfo_o(rinfo_o)
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
    void SpatialSplit::TaskSplitParallel<Prim>::task_split_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      split->split<false>(threadIndex,threadCount,NULL,alloc,scene,prims,lprims_o,linfos[taskIndex],rprims_o,rinfos[taskIndex]);
    }
    
    template<>
    void SpatialSplit::Split::split<true>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, 
					  PrimRefBlockAlloc<BezierPrim>& alloc, Scene* scene, BezierRefList& prims, 
					  BezierRefList& lprims_o, PrimInfo& linfo_o, 
					  BezierRefList& rprims_o, PrimInfo& rinfo_o) const
    {
      TaskSplitParallel<BezierPrim>(threadIndex,threadCount,scheduler,this,alloc,scene,prims,lprims_o,linfo_o,rprims_o,rinfo_o);
    }

    template<>
    void SpatialSplit::Split::split<true>(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, 
					  PrimRefBlockAlloc<PrimRef>& alloc, Scene* scene, TriRefList& prims, 
					  TriRefList& lprims_o, PrimInfo& linfo_o, 
					  TriRefList& rprims_o, PrimInfo& rinfo_o) const
    {
      TaskSplitParallel<PrimRef>(threadIndex,threadCount,scheduler,this,alloc,scene,prims,lprims_o,linfo_o,rprims_o,rinfo_o);
    }
  }
}

  
