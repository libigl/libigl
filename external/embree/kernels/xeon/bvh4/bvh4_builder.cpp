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

#include "bvh4.h"
#include "bvh4_builder.h"
#include "bvh4_refit.h"
#include "bvh4_rotate.h"
#include "bvh4_statistics.h"

#include "builders/heuristics.h"
#include "builders/splitter_fallback.h"

#define ROTATE_TREE 1

#include "common/scene_triangle_mesh.h"

namespace embree
{
  template<typename Heuristic>
  void BVH4Builder<Heuristic>::build(size_t threadIndex, size_t threadCount) 
  {
#if 0
    bvh->clear();
    if (source->isEmpty()) 
      return;
#else
    size_t numPrimitives = source->size();
    bvh->init(numPrimitives);
    if (source->isEmpty()) 
      return;
#endif

    if (g_verbose >= 2) 
      std::cout << "building BVH4<" << bvh->primTy.name << "> with " << Heuristic::name() << " SAH builder ... " << std::flush;

    double t0 = 0.0, t1 = 0.0f;
    if (g_verbose >= 2 || g_benchmark)
      t0 = getSeconds();
    
    /* first generate primrefs */
    new (&initStage) PrimRefGenNormal(threadIndex,threadCount,source,&alloc);
    bvh->numPrimitives = initStage.numPrimitives;
    if (primTy.needVertices) bvh->numVertices = initStage.numVertices;
    else                     bvh->numVertices = 0;

    /* now build BVH */
    TaskScheduler::executeTask(threadIndex,threadCount,_buildFunction,this,"BVH4Builder::build");

    /* finish build */
#if ROTATE_TREE
    for (int i=0; i<5; i++) 
      BVH4Rotate::rotate(bvh,bvh->root);
#endif
    bvh->clearBarrier(bvh->root);
    bvh->bounds = initStage.pinfo.geomBounds;
    initStage.pinfo.clear();

    /* free all temporary blocks */
    Alloc::global.clear();

    if (g_verbose >= 2 || g_benchmark) 
      t1 = getSeconds();
    
    if (g_verbose >= 2) {
      std::cout << "[DONE]" << std::endl;
      std::cout << "  dt = " << 1000.0f*(t1-t0) << "ms, perf = " << 1E-6*double(source->size())/(t1-t0) << " Mprim/s" << std::endl;
      std::cout << BVH4Statistics(bvh).str();
    }

    if (g_benchmark) {
      BVH4Statistics stat(bvh);
      std::cout << "BENCHMARK_BUILD " << 1000.0f*(t1-t0) << " " << 1E-6*double(source->size())/(t1-t0) << " " << stat.bytesUsed() << std::endl;
    }
  }

  template<typename Heuristic>
  BVH4Builder<Heuristic>::BVH4Builder (BVH4* bvh, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize)
    : source(source), geometry(geometry), primTy(bvh->primTy), minLeafSize(minLeafSize), maxLeafSize(maxLeafSize), bvh(bvh),
      taskQueue(Heuristic::depthFirst ? TaskScheduler::GLOBAL_BACK : TaskScheduler::GLOBAL_FRONT)
  {
    size_t maxLeafPrims = BVH4::maxLeafBlocks*primTy.blockSize;
    if (maxLeafPrims < this->maxLeafSize) 
      this->maxLeafSize = maxLeafPrims;
  }
  
  template<typename Heuristic>
  void BVH4Builder<Heuristic>::buildFunction(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event) {
    recurse(threadIndex,threadCount,event,bvh->root,1,initStage.prims,initStage.pinfo,initStage.split);
  }
  
  template<typename Heuristic>
  typename BVH4Builder<Heuristic>::NodeRef BVH4Builder<Heuristic>::createLeaf(size_t threadIndex, atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo)
  {
    /* allocate leaf node */
    size_t blocks = primTy.blocks(pinfo.size());
    char* leaf = bvh->allocPrimitiveBlocks(threadIndex,blocks);
    assert(blocks <= (size_t)BVH4::maxLeafBlocks);

    /* insert all triangles */
    atomic_set<PrimRefBlock>::block_iterator_unsafe iter(prims);
    for (size_t i=0; i<blocks; i++) {
      primTy.pack(leaf+i*primTy.bytes,iter,geometry);
    }
    assert(!iter);
    
    /* free all primitive blocks */
    while (atomic_set<PrimRefBlock>::item* block = prims.take())
      alloc.free(threadIndex,block);
        
    return bvh->encodeLeaf(leaf,blocks);
  }

  template<typename Heuristic>
  typename BVH4Builder<Heuristic>::NodeRef BVH4Builder<Heuristic>::createLargeLeaf(size_t threadIndex, atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, size_t depth)
  {
#if defined(_DEBUG)
    if (depth >= BVH4::maxBuildDepthLeaf) 
      throw std::runtime_error("ERROR: Loosing primitives during build.");
#endif
    
    /* create leaf for few primitives */
    if (pinfo.size() <= maxLeafSize)
      return createLeaf(threadIndex,prims,pinfo);

    /* first level */
    atomic_set<PrimRefBlock> prims0, prims1;
    PrimInfo                 cinfo0, cinfo1;
    FallBackSplitter<Heuristic>::split(threadIndex,&alloc,source,prims,pinfo,prims0,cinfo0,prims1,cinfo1);

    /* second level */
    atomic_set<PrimRefBlock> cprims[4];
    PrimInfo                 cinfo[4];
    FallBackSplitter<Heuristic>::split(threadIndex,&alloc,source,prims0,cinfo0,cprims[0],cinfo[0],cprims[1],cinfo[1]);
    FallBackSplitter<Heuristic>::split(threadIndex,&alloc,source,prims1,cinfo1,cprims[2],cinfo[2],cprims[3],cinfo[3]);

    /*! create an inner node */
    Node* node = bvh->allocNode(threadIndex);
    for (size_t i=0; i<4; i++) node->set(i,cinfo[i].geomBounds,createLargeLeaf(threadIndex,cprims[i],cinfo[i],depth+1));
    return bvh->encodeNode(node);
  }  
  
  template<typename Heuristic>
  void BVH4Builder<Heuristic>::recurse(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event, 
                                       NodeRef& node, size_t depth, atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
  {
    /* use full single threaded build for small jobs */
    if (pinfo.size() < 4*1024) 
      new BuildTask(threadIndex,threadCount,event,this,node,depth,prims,pinfo,split);
    
    /* use single threaded split for medium size jobs  */
    else if (pinfo.size() < 256*1024)
      new SplitTask(threadIndex,threadCount,event,this,node,depth,prims,pinfo,split);
      
    /* use parallel splitter for big jobs */
    else new ParallelSplitTask(threadIndex,threadCount,event,this,node,depth,prims,pinfo,split);
  }

  template<typename Heuristic>
  BVH4Builder<Heuristic>::BuildTask::BuildTask(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event, 
                                               BVH4Builder* parent, NodeRef& node, size_t depth, 
                                               atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
    : threadIndex(threadIndex), parent(parent), dst(node), depth(depth), prims(prims), pinfo(pinfo), split(split)
  {
    new (&task) TaskScheduler::Task(event,_run,this,"build::full");
    TaskScheduler::addTask(threadIndex,parent->taskQueue,&task);
  }
  
  template<typename Heuristic>
  void BVH4Builder<Heuristic>::BuildTask::run(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
  {
    this->threadIndex = threadIndex;
    dst = recurse(depth,prims,pinfo,split);
#if ROTATE_TREE
    for (int i=0; i<5; i++) 
      BVH4Rotate::rotate(parent->bvh,dst); 
#endif
    dst.setBarrier();
    delete this;
  }

  template<typename Heuristic>
  typename BVH4Builder<Heuristic>::NodeRef BVH4Builder<Heuristic>::BuildTask::recurse(size_t depth, atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
  {
    /*! compute leaf and split cost */
    const float leafSAH  = parent->primTy.intCost*pinfo.sah();
    const float splitSAH = BVH4::travCost*halfArea(pinfo.geomBounds)+parent->primTy.intCost*split.sah();
    assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(prims).size() == pinfo.size());
    assert(pinfo.size() == 0 || leafSAH >= 0 && splitSAH >= 0);
    
    /*! create a leaf node when threshold reached or SAH tells us to stop */
    if (pinfo.size() <= parent->minLeafSize || depth > BVH4::maxBuildDepth || (pinfo.size() <= parent->maxLeafSize && leafSAH <= splitSAH)) {
      return parent->createLargeLeaf(threadIndex,prims,pinfo,depth+1);
    }
    
    /*! initialize child list */
    atomic_set<PrimRefBlock> cprims[BVH4::N]; cprims[0] = prims;
    PrimInfo                 cinfo [BVH4::N]; cinfo [0] = pinfo;
    Split                    csplit[BVH4::N]; csplit[0] = split;
    size_t numChildren = 1;
    
    /*! split until node is full or SAH tells us to stop */
    do {
      
      /*! find best child to split */
      float bestSAH = 0; 
      ssize_t bestChild = -1;
      for (size_t i=0; i<numChildren; i++) 
      {
        float dSAH = csplit[i].sah()-cinfo[i].sah();
        if (cinfo[i].size() <= parent->minLeafSize) continue; 
        if (cinfo[i].size() > parent->maxLeafSize) dSAH = min(0.0f,dSAH); //< force split for large jobs
        if (dSAH <= bestSAH) { bestChild = i; bestSAH = dSAH; }
      }
      if (bestChild == -1) break;
      
      /*! perform best found split and find new splits */
      SplitterNormal splitter(threadIndex,&parent->alloc,parent->source,cprims[bestChild],cinfo[bestChild],csplit[bestChild]);
      cprims[bestChild  ] = splitter.lprims; cinfo[bestChild  ] = splitter.linfo; csplit[bestChild  ] = splitter.lsplit;
      cprims[numChildren] = splitter.rprims; cinfo[numChildren] = splitter.rinfo; csplit[numChildren] = splitter.rsplit;
      numChildren++;
      
    } while (numChildren < BVH4::N);
    
    /*! create an inner node */
    Node* node = parent->bvh->allocNode(threadIndex);
    for (size_t i=0; i<numChildren; i++) node->set(i,cinfo[i].geomBounds,recurse(depth+1,cprims[i],cinfo[i],csplit[i]));
    return parent->bvh->encodeNode(node);
  }
  
  template<typename Heuristic>
  BVH4Builder<Heuristic>::SplitTask::SplitTask(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event, 
                                               BVH4Builder* parent, NodeRef& node, size_t depth, 
                                               atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
    : parent(parent), dst(node), depth(depth), prims(prims), pinfo(pinfo), split(split)
  {
    new (&task) TaskScheduler::Task(event,_recurse,this,"build::split");
    TaskScheduler::addTask(threadIndex,parent->taskQueue,&task);
  }
  
  template<typename Heuristic>
  void BVH4Builder<Heuristic>::SplitTask::recurse(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
  {
    /*! compute leaf and split cost */
    const float leafSAH  = parent->primTy.intCost*pinfo.sah();
    const float splitSAH = BVH4::travCost*halfArea(pinfo.geomBounds)+parent->primTy.intCost*split.sah();
    assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(prims).size() == pinfo.size());
    assert(pinfo.size() == 0 || leafSAH >= 0 && splitSAH >= 0);
    
    /*! create a leaf node when threshold reached or SAH tells us to stop */
    if (pinfo.size() <= parent->minLeafSize || depth > BVH4::maxBuildDepth || (pinfo.size() <= parent->maxLeafSize && leafSAH <= splitSAH)) {
      dst = parent->createLargeLeaf(threadIndex,prims,pinfo,depth+1); delete this; return;
    }
    
    /*! initialize child list */
    atomic_set<PrimRefBlock> cprims[BVH4::N]; 
    PrimInfo                 cinfo [BVH4::N]; 
    Split                    csplit[BVH4::N]; 
    cprims[0] = prims;
    cinfo [0] = pinfo;
    csplit[0] = split;
    size_t numChildren = 1;
    
    /*! split until node is full or SAH tells us to stop */
    do {
      
      /*! find best child to split */
      float bestSAH = 0; 
      ssize_t bestChild = -1;
      for (size_t i=0; i<numChildren; i++) 
      {
        float dSAH = csplit[i].sah()-cinfo[i].sah();
        if (cinfo[i].size() <= parent->minLeafSize) continue; 
        if (cinfo[i].size() > parent->maxLeafSize) dSAH = min(0.0f,dSAH); //< force split for large jobs
        if (dSAH <= bestSAH) { bestChild = i; bestSAH = dSAH; }
      }
      if (bestChild == -1) break;
      
      /*! perform best found split and find new splits */
      SplitterNormal splitter(threadIndex,&parent->alloc,parent->source,cprims[bestChild],cinfo[bestChild],csplit[bestChild]);
      cprims[bestChild  ] = splitter.lprims; cinfo[bestChild  ] = splitter.linfo; csplit[bestChild  ] = splitter.lsplit;
      cprims[numChildren] = splitter.rprims; cinfo[numChildren] = splitter.rinfo; csplit[numChildren] = splitter.rsplit;
      numChildren++;
      
    } while (numChildren < BVH4::N);
    
    /*! create an inner node */
    Node* node = parent->bvh->allocNode(threadIndex);
    dst = parent->bvh->encodeNode(node);
    for (size_t i=0; i<numChildren; i++) {
      node->set(i,cinfo[i].geomBounds,0);
      parent->recurse(threadIndex,threadCount,event,node->child(i),depth+1,cprims[i],cinfo[i],csplit[i]);
    }
    delete this;
  }
  
  template<typename Heuristic>
  BVH4Builder<Heuristic>::ParallelSplitTask::ParallelSplitTask (size_t threadIndex, size_t threadCount, TaskScheduler::Event* event, 
                                                                BVH4Builder* parent, NodeRef& node, size_t depth, 
                                                                atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
    : parent(parent), dst(node), depth(depth), numChildren(1), bestChild(0)
  {
    /*! compute leaf and split cost */
    const float leafSAH  = parent->primTy.intCost*pinfo.sah();
    const float splitSAH = BVH4::travCost*halfArea(pinfo.geomBounds)+parent->primTy.intCost*split.sah();
    assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(prims).size() == pinfo.size());
    assert(pinfo.size() == 0 || leafSAH >= 0 && splitSAH >= 0);
    
    /*! create a leaf node when threshold reached or SAH tells us to stop */
    if (pinfo.size() <= parent->minLeafSize || depth > BVH4::maxBuildDepth || (pinfo.size() <= parent->maxLeafSize && leafSAH <= splitSAH)) {
      dst = parent->createLargeLeaf(threadIndex,prims,pinfo,depth+1); delete this; return;
    }
    
    /*! initialize child list */
    cprims[0] = prims; cinfo [0] = pinfo; csplit[0] = split; 
    
    /*! perform first split */
    new (&splitter) MultiThreadedSplitterNormal(threadIndex,threadCount,event,
                                                &parent->alloc,parent->source,
                                                cprims[bestChild],cinfo[bestChild],csplit[bestChild],
                                                _loop,this);
  }
  
  template<typename Heuristic>
  void BVH4Builder<Heuristic>::ParallelSplitTask::loop(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
  {
    /*! copy new children into work array */
    cprims[bestChild  ] = splitter.lprims; cinfo[bestChild  ] = splitter.linfo; csplit[bestChild  ] = splitter.lsplit;
    cprims[numChildren] = splitter.rprims; cinfo[numChildren] = splitter.rinfo; csplit[numChildren] = splitter.rsplit;
    numChildren++;
    
    /*! find best child to split */
    if (numChildren < BVH4::N) 
    {
      float bestSAH = 0; 
      bestChild = -1;
      for (size_t i=0; i<numChildren; i++) 
      {
        float dSAH = csplit[i].sah()-cinfo[i].sah();
        if (cinfo[i].size() <= parent->minLeafSize) continue; 
        if (cinfo[i].size() > parent->maxLeafSize) dSAH = min(0.0f,dSAH); //< force split for large jobs
        if (dSAH <= bestSAH) { bestChild = i; bestSAH = dSAH; }
      }
      
      /*! perform next split */
      if (bestChild != -1) 
      {
        new (&splitter) MultiThreadedSplitterNormal(threadIndex,threadCount,event,
                                                    &parent->alloc,parent->source,
                                                    cprims[bestChild],cinfo[bestChild],csplit[bestChild],
                                                    _loop,this);
        return;
      }
    }
    
    /*! create an inner node */
    Node* node = parent->bvh->allocNode(threadIndex);
    dst = parent->bvh->encodeNode(node);
    for (size_t i=0; i<numChildren; i++) {
      node->set(i,cinfo[i].geomBounds,0);
      parent->recurse(threadIndex,threadCount,event,node->child(i),depth+1,cprims[i],cinfo[i],csplit[i]);
    }
    delete this;
  }
  
  Builder* BVH4BuilderObjectSplit1 (void* accel, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize) {
    return new BVH4Builder<HeuristicBinning<0> >((BVH4*)accel,source,geometry,minLeafSize,maxLeafSize);
  }

  Builder* BVH4BuilderObjectSplit4 (void* accel, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize) {
    return new BVH4Builder<HeuristicBinning<2> >((BVH4*)accel,source,geometry,minLeafSize,maxLeafSize);
  }

  Builder* BVH4BuilderObjectSplit8 (void* accel, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize) {
    return new BVH4Builder<HeuristicBinning<3> >((BVH4*)accel,source,geometry,minLeafSize,maxLeafSize);
  }

  Builder* BVH4BuilderSpatialSplit1 (void* accel, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize) {
    return new BVH4Builder<HeuristicSpatial<0> >((BVH4*)accel,source,geometry,minLeafSize,maxLeafSize);
  }

  Builder* BVH4BuilderSpatialSplit4 (void* accel, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize) {
    return new BVH4Builder<HeuristicSpatial<2> >((BVH4*)accel,source,geometry,minLeafSize,maxLeafSize);
  }

  Builder* BVH4BuilderSpatialSplit8 (void* accel, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize) {
    return new BVH4Builder<HeuristicSpatial<3> >((BVH4*)accel,source,geometry,minLeafSize,maxLeafSize);
  }
}
