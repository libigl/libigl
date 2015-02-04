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

#include "bvh4.h"
#include "bvh4_builder_hair.h"
#include "bvh4_statistics.h"
#include "common/scene_bezier_curves.h"
#include "../builders/bezierrefgen.h"
#include <algorithm>
#include "geometry/bezier1v.h"
#include "geometry/bezier1i.h"

namespace embree
{
  namespace isa
  {
    template<> BVH4BuilderHairT<Bezier1v>::BVH4BuilderHairT (BVH4* bvh, Scene* scene, size_t mode) : BVH4BuilderHair(bvh,scene,mode) {}
    template<> BVH4BuilderHairT<Bezier1i>::BVH4BuilderHairT (BVH4* bvh, Scene* scene, size_t mode) : BVH4BuilderHair(bvh,scene,mode) {}

    BVH4BuilderHair::BVH4BuilderHair (BVH4* bvh, Scene* scene, size_t mode)
      : scene(scene), minLeafSize(1), maxLeafSize(inf), enableSpatialSplits(mode & MODE_HIGH_QUALITY), listMode(mode & LIST_MODE_BITS), bvh(bvh), scheduler(&scene->lockstep_scheduler), remainingReplications(0)
    {
      if (BVH4::maxLeafBlocks < this->maxLeafSize) 
	this->maxLeafSize = BVH4::maxLeafBlocks;
    }
    
    void BVH4BuilderHair::build(size_t threadIndex, size_t threadCount) 
    {
      /* fast path for empty BVH */
      size_t numPrimitives = scene->numBezierCurves;
      size_t maxPrimitives = max(numPrimitives,(size_t)(g_hair_builder_replication_factor*numPrimitives));
      bvh->init(sizeof(BVH4::UnalignedNode),maxPrimitives,threadCount);
      if (numPrimitives == 0) return;
      numGeneratedPrims = 0;
      
      //while (true)
      {
	double t0 = 0.0;
	if (g_verbose >= 2) 
	{
	  std::cout << "building BVH4<" + bvh->primTy.name + "> using " << TOSTRING(isa) << "::BVH4BuilderHair ..." << std::flush;
	  t0 = getSeconds();
	}
	
	size_t N = 0;
	float r = 0;
	
	/* create initial curve list */
	size_t numVertices = 0;
	BezierRefGen gen(threadIndex,threadCount,scheduler,&alloc,scene);
	PrimInfo pinfo = gen.pinfo;
	BezierRefList prims = gen.prims;
	
	bvh->numPrimitives = scene->numBezierCurves;
	bvh->numVertices = 0;
	if (&bvh->primTy == &SceneBezier1i::type) bvh->numVertices = numVertices;
	
	/* start recursive build */
	remainingReplications = maxPrimitives-numPrimitives;
	bvh->bounds = pinfo.geomBounds;

	/* return if all geometry got filtered out */
	if (gen.pinfo.size() == 0)
	  return;

        Allocator nodeAlloc(&bvh->alloc);
        Allocator leafAlloc(&bvh->alloc);
	
#if 0
	const Split split = find_split(threadIndex,threadCount,prims,pinfo,pinfo.geomBounds);
	BuildTask task(&bvh->root,0,prims,pinfo,pinfo.geomBounds,split); recurseTask(threadIndex,nodeAlloc,leafAlloc,task);
	_mm_sfence(); // make written leaves globally visible
#else
	const Split split = find_split<true>(threadIndex,threadCount,prims,pinfo,pinfo.geomBounds,pinfo);
	BuildTask task(&bvh->root,0,prims,pinfo,pinfo.geomBounds,pinfo,split);
	numActiveTasks = 1;
	tasks.push_back(task);
	std::push_heap(tasks.begin(),tasks.end());
	
#if 1
	while (tasks.front().pinfo.size() > 200000)
	{
	  BuildTask task = tasks.front();
	  std::pop_heap(tasks.begin(),tasks.end());
	  tasks.pop_back();
	  
	  size_t numChildren;
	  BuildTask ctasks[BVH4::N];
	  processTask<true>(threadIndex,threadCount,nodeAlloc,leafAlloc,task,ctasks,numChildren);
	  
	  for (size_t i=0; i<numChildren; i++) {
	    atomic_add(&numActiveTasks,+1);
	    tasks.push_back(ctasks[i]);
	    std::push_heap(tasks.begin(),tasks.end());
	  }
	  atomic_add(&numActiveTasks,-1);
	}
	_mm_sfence(); // make written leaves globally visible
#endif
	
	scheduler->dispatchTask(threadIndex,threadCount,_task_build_parallel,this,threadCount,"BVH4Builder::build_parallel");

        tasks.clear();
#endif
	
	if (g_verbose >= 2) {
	  double t1 = getSeconds();
	  std::cout << " [DONE]" << std::endl;
	  std::cout << "  dt = " << 1000.0f*(t1-t0) << "ms, perf = " << 1E-6*double(numPrimitives)/(t1-t0) << " Mprim/s" << std::endl;
	  std::cout << BVH4Statistics(bvh).str();
	}
      }
    }

    template<typename Primitive>
    BVH4::NodeRef BVH4BuilderHairT<Primitive>::createLeaf(size_t threadIndex, Allocator& nodeAlloc, Allocator& leafAlloc, 
                                                          size_t depth, BezierRefList& prims, const PrimInfo& pinfo)
    {
      size_t N = pinfo.size();
      
      if (N > (size_t)BVH4::maxLeafBlocks) {
	//std::cout << "WARNING: Loosing " << N-BVH4::maxLeafBlocks << " primitives during build!" << std::endl;
	std::cout << "!" << std::flush;
	N = (size_t)BVH4::maxLeafBlocks;
      }
      if (g_verbose >= 1) {
        size_t numGeneratedPrimsOld = atomic_add(&numGeneratedPrims,N); 
        if (numGeneratedPrimsOld%10000 > (numGeneratedPrimsOld+N)%10000) std::cout << "." << std::flush; 
      }
      //assert(N <= (size_t)BVH4::maxLeafBlocks);
     
      Primitive* leaf = (Primitive*) bvh->allocPrimitiveBlocks(nodeAlloc,N);
      BezierRefList::block_iterator_unsafe iter(prims);
      for (size_t i=0; i<N; i++) leaf[i].fill(iter,scene,listMode);
      assert(!iter);
      
      /* free all primitive blocks */
      while (BezierRefList::item* block = prims.take())
	alloc.free(threadIndex,block);
      
      return bvh->encodeLeaf((char*)leaf,listMode ? listMode : N);
    }

    BVH4::NodeRef BVH4BuilderHair::createLargeLeaf(size_t threadIndex, Allocator& nodeAlloc, Allocator& leafAlloc, BezierRefList& prims, const PrimInfo& pinfo, size_t depth)
    {
#if defined(_DEBUG)
      if (depth >= BVH4::maxBuildDepthLeaf) 
	THROW_RUNTIME_ERROR("ERROR: Loosing primitives during build.");
#endif
      
      /* create leaf for few primitives */
      if (pinfo.size() <= BVH4::maxLeafBlocks)
	return createLeaf(threadIndex,nodeAlloc,leafAlloc,depth,prims,pinfo);
      
      /* first level */
      BezierRefList prims0, prims1;
      PrimInfo   cinfo0, cinfo1;
      FallBackSplit::find(threadIndex,alloc,prims,prims0,cinfo0,prims1,cinfo1);
      
      /* second level */
      BezierRefList cprims[4];
      PrimInfo   cinfo[4];
      FallBackSplit::find(threadIndex,alloc,prims0,cprims[0],cinfo[0],cprims[1],cinfo[1]);
      FallBackSplit::find(threadIndex,alloc,prims1,cprims[2],cinfo[2],cprims[3],cinfo[3]);
      
      /*! create an inner node */
      BVH4::Node* node = bvh->allocNode(nodeAlloc);
      for (size_t i=0; i<4; i++) {
        node->set(i,cinfo[i].geomBounds);
        node->set(i,createLargeLeaf(threadIndex,nodeAlloc,leafAlloc,cprims[i],cinfo[i],depth+1));
      }
      BVH4::compact(node); // moves empty nodes to the end
      return bvh->encodeNode(node);
    }  

    template<bool Parallel>
    Split BVH4BuilderHair::find_split(size_t threadIndex, size_t threadCount, BezierRefList& prims, const PrimInfo& pinfo, const NAABBox3fa& bounds, const PrimInfo& sinfo)
    {
      /* variable to track the SAH of the best splitting approach */
      float bestSAH = inf;
      const float leafSAH = BVH4::intCost*float(pinfo.size())*halfArea(bounds.bounds);
      
      /* perform standard binning in aligned space */
      ObjectPartition::Split alignedObjectSplit;
      float alignedObjectSAH = inf;
      alignedObjectSplit = ObjectPartition::find<Parallel>(threadIndex,threadCount,scheduler,prims,pinfo,0); // FIXME: hardcoded 0
      alignedObjectSAH = BVH4::travCostAligned*halfArea(bounds.bounds) + BVH4::intCost*alignedObjectSplit.splitSAH();
      bestSAH = min(alignedObjectSAH,bestSAH);
      
      /* perform spatial split in aligned space */
      SpatialSplit::Split alignedSpatialSplit;
      float alignedSpatialSAH = inf;
      if (enableSpatialSplits && remainingReplications > 0) {
	alignedSpatialSplit = SpatialSplit::find<Parallel>(threadIndex,threadCount,scheduler,scene,prims,pinfo,0); // FIXME: hardcoded 0
	alignedSpatialSAH = BVH4::travCostAligned*halfArea(bounds.bounds) + BVH4::intCost*alignedSpatialSplit.splitSAH();
	bestSAH = min(alignedSpatialSAH,bestSAH);
      }
      
      /* perform standard binning in unaligned space */
      ObjectPartitionUnaligned::Split unalignedObjectSplit;
      float unalignedObjectSAH = inf;
      if (alignedObjectSAH > 0.7f*leafSAH) {
	if (sinfo.size()) 
	  unalignedObjectSplit = ObjectPartitionUnaligned::find<Parallel>(threadIndex,threadCount,scheduler,prims,bounds.space,sinfo);
	else {
	  const LinearSpace3fa space = ObjectPartitionUnaligned::computeAlignedSpace(threadIndex,threadCount,scheduler,prims); 
	  const PrimInfo       sinfo = ObjectPartitionUnaligned::computePrimInfo    <Parallel>(threadIndex,threadCount,scheduler,prims,space);
	  unalignedObjectSplit = ObjectPartitionUnaligned::find<Parallel>(threadIndex,threadCount,scheduler,prims,space,sinfo);
	}    	
	unalignedObjectSAH = BVH4::travCostUnaligned*halfArea(bounds.bounds) + BVH4::intCost*unalignedObjectSplit.splitSAH();
	bestSAH = min(unalignedObjectSAH,bestSAH);
      }
      
      /* perform splitting into two strands */
      StrandSplit::Split strandSplit;
      float strandSAH = inf;
      if (alignedObjectSAH > 0.6f*leafSAH) {
	strandSplit = StrandSplit::find<Parallel>(threadIndex,threadCount,scheduler,prims);
	strandSAH = BVH4::travCostUnaligned*halfArea(bounds.bounds) + BVH4::intCost*strandSplit.splitSAH();
	bestSAH = min(strandSAH,bestSAH);
      }

      
      /* return best split */
      if      (bestSAH == float(inf)        ) return Split();
      else if (bestSAH == alignedObjectSAH  ) return Split(alignedObjectSplit);
      else if (bestSAH == alignedSpatialSAH ) return Split(alignedSpatialSplit);
      else if (bestSAH == unalignedObjectSAH) return Split(unalignedObjectSplit);
      else if (bestSAH == strandSAH         ) return Split(strandSplit);
      else THROW_RUNTIME_ERROR("bvh4hair_builder: internal error");
    }
    
    template<bool Parallel>
    __forceinline void BVH4BuilderHair::processTask(size_t threadIndex, size_t threadCount, Allocator& nodeAlloc, Allocator& leafAlloc,
                                                    BuildTask& task, BuildTask task_o[BVH4::N], size_t& numTasks_o)
    {
      /* create enforced leaf */
      const float leafSAH  = BVH4::intCost*task.pinfo.leafSAH();
      const float splitSAH = BVH4::travCostUnaligned*halfArea(task.bounds.bounds)+BVH4::intCost*task.split.splitSAH();
      
      if (task.pinfo.size() <= minLeafSize || task.depth >= BVH4::maxBuildDepth || (task.pinfo.size() <= maxLeafSize && leafSAH <= splitSAH)) {
	*task.dst = createLargeLeaf(threadIndex,nodeAlloc,leafAlloc,task.prims,task.pinfo,task.depth);
	numTasks_o = 0;
	return;
      }
      
      /*! initialize child list */
      bool isAligned = true;
      PrimInfo cpinfo     [BVH4::N]; cpinfo [0] = task.pinfo; 
      PrimInfo csinfo     [BVH4::N]; csinfo [0] = task.sinfo; 
      NAABBox3fa cbounds  [BVH4::N]; cbounds[0] = task.bounds;
      BezierRefList cprims[BVH4::N]; cprims [0] = task.prims;
      Split csplit        [BVH4::N]; csplit [0] = task.split;        
      size_t numChildren = 1;
      
      /*! split until node is full or SAH tells us to stop */
      do {
	
	/*! find best child to split */
	float bestSAH = 0; 
	ssize_t bestChild = -1;
	for (size_t i=0; i<numChildren; i++) 
	{
	  float dSAH = csplit[i].splitSAH()-cpinfo[i].leafSAH();
	  if (cpinfo[i].size() <= minLeafSize) continue; 
	  if (cpinfo[i].size() > maxLeafSize) dSAH = min(0.0f,dSAH); //< force split for large jobs
	  if (dSAH <= bestSAH) { bestChild = i; bestSAH = dSAH; }
	}
	if (bestChild == -1) break;
	
	/*! split selected child */
	PrimInfo linfo(empty), rinfo(empty);
	BezierRefList lprims, rprims;
	csplit[bestChild].split<Parallel>(threadIndex,threadCount,scheduler,alloc,scene,cprims[bestChild],lprims,linfo,rprims,rinfo);
	const ssize_t replications = linfo.size()+rinfo.size()-cpinfo[bestChild].size(); assert(replications >= 0);
	isAligned &= csplit[bestChild].isAligned;
	LinearSpace3fa lspace,rspace;
	PrimInfo lsinfo(empty),rsinfo(empty);
	if (!isAligned) {
	  lspace = ObjectPartitionUnaligned::computeAlignedSpace(threadIndex,threadCount,scheduler,lprims); 
	  rspace = ObjectPartitionUnaligned::computeAlignedSpace(threadIndex,threadCount,scheduler,rprims); 
	  lsinfo = ObjectPartitionUnaligned::computePrimInfo    <Parallel>(threadIndex,threadCount,scheduler,lprims,lspace);
	  rsinfo = ObjectPartitionUnaligned::computePrimInfo    <Parallel>(threadIndex,threadCount,scheduler,rprims,rspace);
	}
	const NAABBox3fa lbounds = isAligned ? linfo.geomBounds : NAABBox3fa(lspace,lsinfo.geomBounds); 
	const NAABBox3fa rbounds = isAligned ? rinfo.geomBounds : NAABBox3fa(rspace,rsinfo.geomBounds); 
	const Split lsplit = find_split<Parallel>(threadIndex,threadCount,lprims,linfo,lbounds,lsinfo);
	const Split rsplit = find_split<Parallel>(threadIndex,threadCount,rprims,rinfo,rbounds,rsinfo);
	cprims[numChildren] = rprims; cpinfo[numChildren] = rinfo; cbounds[numChildren] = rbounds; csinfo[numChildren] = lsinfo; csplit[numChildren] = rsplit;
	cprims[bestChild  ] = lprims; cpinfo[bestChild  ] = linfo; cbounds[bestChild  ] = lbounds; csinfo[bestChild  ] = rsinfo; csplit[bestChild  ] = lsplit;
	if (replications) atomic_add(&remainingReplications,-replications); 
	numChildren++;
	
      } while (numChildren < BVH4::N);
      
      /* create aligned node */
      if (isAligned) 
      {
	BVH4::Node* node = bvh->allocNode(nodeAlloc);
	for (size_t i=0; i<numChildren; i++) {
	  node->set(i,cpinfo[i].geomBounds);
	  new (&task_o[i]) BuildTask(&node->child(i),task.depth+1,cprims[i],cpinfo[i],cbounds[i],csinfo[i],csplit[i]);
	}
	numTasks_o = numChildren;
	*task.dst = bvh->encodeNode(node);
      }
      
      /* create unaligned node */
      else {
	BVH4::UnalignedNode* node = bvh->allocUnalignedNode(nodeAlloc);
	for (size_t i=0; i<numChildren; i++) {
	  node->set(i,cbounds[i]);
	  new (&task_o[i]) BuildTask(&node->child(i),task.depth+1,cprims[i],cpinfo[i],cbounds[i],csinfo[i],csplit[i]);
	}
	numTasks_o = numChildren;
	*task.dst = bvh->encodeNode(node);
      }
    }
    
    void BVH4BuilderHair::recurseTask(size_t threadIndex, size_t threadCount, Allocator& nodeAlloc, Allocator& leafAlloc, BuildTask& task)
    {
      size_t numChildren;
      BuildTask tasks[BVH4::N];
      processTask<false>(threadIndex,threadCount,nodeAlloc,leafAlloc,task,tasks,numChildren);
      for (size_t i=0; i<numChildren; i++) 
	recurseTask(threadIndex,threadCount,nodeAlloc,leafAlloc,tasks[i]);
    }
    
    void BVH4BuilderHair::task_build_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      Allocator nodeAlloc(&bvh->alloc);
      Allocator leafAlloc(&bvh->alloc);

      while (numActiveTasks) 
      {
	taskMutex.lock();
	if (tasks.size() == 0) {
	  taskMutex.unlock();
	  continue;
	}
	
	/* take next task from heap */
	BuildTask task = tasks.back(); tasks.pop_back();
	taskMutex.unlock();
	
	/* recursively finish task */
	if (task.pinfo.size() < 1024) {
	  atomic_add(&numActiveTasks,-1);
	  recurseTask(threadIndex,threadCount,nodeAlloc,leafAlloc,task);
	}
	
	/* execute task and add child tasks */
	else 
	{
	  size_t numChildren;
	  BuildTask ctasks[BVH4::N];
	  processTask<false>(threadIndex,threadCount,nodeAlloc,leafAlloc,task,ctasks,numChildren);
	  taskMutex.lock();
	  for (size_t i=0; i<numChildren; i++) {
	    atomic_add(&numActiveTasks,+1);
	    tasks.push_back(ctasks[i]);
	  }
	  atomic_add(&numActiveTasks,-1);
	  taskMutex.unlock();
	}
      }
      _mm_sfence(); // make written leaves globally visible
    }

    /*! entry functions for the builder */
    Builder* BVH4Bezier1vBuilder_OBB (void* bvh, Scene* scene, size_t mode) { return new class BVH4BuilderHairT<Bezier1v> ((BVH4*)bvh,scene,mode); }
    Builder* BVH4Bezier1iBuilder_OBB (void* bvh, Scene* scene, size_t mode) { return new class BVH4BuilderHairT<Bezier1i> ((BVH4*)bvh,scene,mode); }
  }
}
