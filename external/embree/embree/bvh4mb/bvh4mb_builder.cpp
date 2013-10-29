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

#include "bvh4mb_builder.h"
#include "../builders/heuristics.h"
#include "../geometry/triangles.h"

namespace embree
{
  template<typename Heuristic>
  BVH4MBBuilder<Heuristic>::BVH4MBBuilder(TaskScheduler::Event* event, BVH4MB* bvh, const size_t minLeafSize, const size_t maxLeafSize)
    : geom(bvh->geom), trity(bvh->trity),
      minLeafSize(minLeafSize), maxLeafSize(maxLeafSize),
      taskQueue(Heuristic::depthFirst ? TaskScheduler::GLOBAL_FRONT : TaskScheduler::GLOBAL_BACK),
      bvh(bvh)
  {
    if (bvh->maxLeafTris < this->maxLeafSize) 
      this->maxLeafSize = bvh->maxLeafTris;

    bvh->clear();
    new (&finishStage) TaskScheduler::EventScheduleTask(event,_finish,this,"BVH4MBBuilder::finish");
    new (&buildStage ) TaskScheduler::EventScheduleTask(&finishStage,_build,this,"BVH4MBBuilder::build");
    new (&initStage  ) PrimRefGenNormal(&buildStage,geom,&alloc);
  }

  template<typename Heuristic>
  BVH4MBBuilder<Heuristic>::~BVH4MBBuilder() {
  }

  template<typename Heuristic>
  void BVH4MBBuilder<Heuristic>::build(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event) {
    recurse(threadIndex,threadCount,event,bvh->root,1,initStage.prims,initStage.pinfo,initStage.split);
  }

  template<typename Heuristic>
  void BVH4MBBuilder<Heuristic>::finish(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event) {
    bvh->refit(geom,bvh->root);
    bvh->bounds = initStage.pinfo.geomBounds;
    delete this;
  }

  template<typename Heuristic>
  void BVH4MBBuilder<Heuristic>::recurse(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event, 
                                         BVH4MB::Base*& node, size_t depth, 
                                         atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
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

  /***********************************************************************************************************************
   *                                         Leaf creation
   **********************************************************************************************************************/

  template<typename Heuristic>
  BVH4MB::Base* BVH4MBBuilder<Heuristic>::createLeaf(size_t threadIndex, atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo)
  {
    /* allocate leaf node */
    size_t blocks = trity.blocks(pinfo.size());
    char* leaf = (char*) bvh->alloc.malloc(threadIndex,blocks*trity.bytes,1 << BVH4MB::alignment);

    /* insert all triangles */
    atomic_set<PrimRefBlock>::block_iterator_unsafe iter(prims);
    for (size_t i=0; i<blocks; i++) trity.pack(leaf+i*trity.bytes,iter,geom);
    assert(!iter);

    /* free all primitive blocks */
    while (atomic_set<PrimRefBlock>::item* block = prims.take())
      alloc.free(threadIndex,block);

    return BVH4MB::Base::encodeLeaf(leaf,blocks);
  }

  /***********************************************************************************************************************
   *                                         Full Recursive Build Task
   **********************************************************************************************************************/

  template<typename Heuristic>
  __forceinline BVH4MBBuilder<Heuristic>::BuildTask::BuildTask(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event, 
                                                               BVH4MBBuilder* parent, BVH4MB::Base*& node, size_t depth, 
                                                               atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
    : threadIndex(threadIndex), parent(parent), dst(node), depth(depth), prims(prims), pinfo(pinfo), split(split)
  {
    new (&task) TaskScheduler::Task(event,_run,this,"build::full");
    TaskScheduler::addTask(threadIndex,parent->taskQueue,&task);
  }

  template<typename Heuristic>
  void BVH4MBBuilder<Heuristic>::BuildTask::run(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
  {
    this->threadIndex = threadIndex;
    dst = recurse(depth,prims,pinfo,split);
    for (int i=0; i<5; i++) parent->bvh->rotate(dst,depth);
    parent->bvh->sort(parent->geom,dst,inf);
    parent->bvh->refit(parent->geom,dst);
    delete this;
  }

  template<typename Heuristic>
  BVH4MB::Base* BVH4MBBuilder<Heuristic>::BuildTask::recurse(size_t depth, atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
  {
    /*! compute leaf and split cost */
    const float leafSAH  = parent->trity.intCost*pinfo.sah();
    const float splitSAH = BVH4MB::travCost*halfArea(pinfo.geomBounds)+parent->trity.intCost*split.sah();
    assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(prims).size() == pinfo.size());
    assert(pinfo.size() == 0 || leafSAH >= 0 && splitSAH >= 0);
    
    /*! create a leaf node when threshold reached or SAH tells us to stop */
    if (pinfo.size() <= parent->minLeafSize || depth > BVH4MB::maxDepth || (pinfo.size() <= parent->maxLeafSize && leafSAH <= splitSAH)) {
      return parent->createLeaf(threadIndex,prims,pinfo);
    }

    /*! initialize child list */
    atomic_set<PrimRefBlock> cprims[4]; cprims[0] = prims;
    PrimInfo                 cinfo [4]; cinfo [0] = pinfo;
    Split                    csplit[4]; csplit[0] = split;
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
      SplitterNormal splitter(threadIndex,&parent->alloc,parent->geom,cprims[bestChild],cinfo[bestChild],csplit[bestChild]);
      cprims[bestChild  ] = splitter.lprims; cinfo[bestChild  ] = splitter.linfo; csplit[bestChild  ] = splitter.lsplit;
      cprims[numChildren] = splitter.rprims; cinfo[numChildren] = splitter.rinfo; csplit[numChildren] = splitter.rsplit;
      numChildren++;
      
    } while (numChildren < 4);

    /*! create an inner node */
    BVH4MB::Node* node = (BVH4MB::Node*) parent->bvh->alloc.malloc(threadIndex,sizeof(BVH4MB::Node),1 << BVH4MB::alignment); node->clear();
    for (size_t i=0; i<numChildren; i++) node->set(i,cinfo[i].geomBounds,recurse(depth+1,cprims[i],cinfo[i],csplit[i]));
    return BVH4MB::Base::encodeNode(node);
  }

  /***********************************************************************************************************************
   *                                      Single Threaded Split Task
   **********************************************************************************************************************/

  template<typename Heuristic>
  __forceinline BVH4MBBuilder<Heuristic>::SplitTask::SplitTask(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event,
                                                               BVH4MBBuilder* parent, BVH4MB::Base*& node, size_t depth, 
                                                               atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
    : parent(parent), dst(node), depth(depth), prims(prims), pinfo(pinfo), split(split)
  {
    new (&task) TaskScheduler::Task(event,_recurse,this,"build::split");
    TaskScheduler::addTask(threadIndex,parent->taskQueue,&task);
  }

  template<typename Heuristic>
  void BVH4MBBuilder<Heuristic>::SplitTask::recurse(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
  {
     /*! compute leaf and split cost */
    const float leafSAH  = parent->trity.intCost*pinfo.sah();
    const float splitSAH = BVH4MB::travCost*halfArea(pinfo.geomBounds)+parent->trity.intCost*split.sah();
    assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(prims).size() == pinfo.size());
    assert(pinfo.size() == 0 || leafSAH >= 0 && splitSAH >= 0);

    /*! create a leaf node when threshold reached or SAH tells us to stop */
    if (pinfo.size() <= parent->minLeafSize || depth > BVH4MB::maxDepth || (pinfo.size() <= parent->maxLeafSize && leafSAH <= splitSAH)) {
      dst = parent->createLeaf(threadIndex,prims,pinfo); delete this; return;
    }

    /*! initialize child list */
    atomic_set<PrimRefBlock> cprims[4]; cprims[0] = prims;
    PrimInfo                 cinfo [4]; cinfo [0] = pinfo;
    Split                    csplit[4]; csplit[0] = split;
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
      SplitterNormal splitter(threadIndex,&parent->alloc,parent->geom,cprims[bestChild],cinfo[bestChild],csplit[bestChild]);
      cprims[bestChild  ] = splitter.lprims; cinfo[bestChild  ] = splitter.linfo; csplit[bestChild  ] = splitter.lsplit;
      cprims[numChildren] = splitter.rprims; cinfo[numChildren] = splitter.rinfo; csplit[numChildren] = splitter.rsplit;
      numChildren++;
      
    } while (numChildren < 4);
     
    /*! create an inner node */
    BVH4MB::Node* node = (BVH4MB::Node*) parent->bvh->alloc.malloc(threadIndex,sizeof(BVH4MB::Node),1 << BVH4MB::alignment); node->clear();
    dst = BVH4MB::Base::encodeNode(node);
    for (size_t i=0; i<numChildren; i++) {
      node->set(i,cinfo[i].geomBounds,NULL);
      parent->recurse(threadIndex,threadCount,event,node->child[i],depth+1,cprims[i],cinfo[i],csplit[i]);
    }
    delete this;
  }
  
  /***********************************************************************************************************************
   *                                           Parallel Split Task
   **********************************************************************************************************************/

  template<typename Heuristic>
  BVH4MBBuilder<Heuristic>::ParallelSplitTask::ParallelSplitTask(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event, 
                                                                 BVH4MBBuilder* parent, BVH4MB::Base*& node, size_t depth, 
                                                                 atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
    : parent(parent), dst(node), depth(depth), numChildren(1), bestChild(0)
  {
    /*! compute leaf and split cost */
    const float leafSAH  = parent->trity.intCost*pinfo.sah();
    const float splitSAH = BVH4MB::travCost*halfArea(pinfo.geomBounds)+parent->trity.intCost*split.sah();
    assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(prims).size() == pinfo.size());
    assert(pinfo.size() == 0 || leafSAH >= 0 && splitSAH >= 0);

    /*! create a leaf node when threshold reached or SAH tells us to stop */
    if (pinfo.size() <= parent->minLeafSize || depth > BVH4MB::maxDepth || (pinfo.size() <= parent->maxLeafSize && leafSAH <= splitSAH)) {
      dst = parent->createLeaf(threadIndex,prims,pinfo); delete this; return;
    }

    /*! initialize child list */
    cprims[0] = prims; cinfo [0] = pinfo; csplit[0] = split; 

    /*! perform first split */
    new (&splitter) MultiThreadedSplitterNormal(threadIndex,threadCount,event,
                                                &parent->alloc,parent->geom,
                                                cprims[bestChild],cinfo[bestChild],csplit[bestChild],
                                                _loop,this);
  }

  template<typename Heuristic>
  void BVH4MBBuilder<Heuristic>::ParallelSplitTask::loop(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
  {
    /*! copy new children into work array */
    cprims[bestChild  ] = splitter.lprims; cinfo[bestChild  ] = splitter.linfo; csplit[bestChild  ] = splitter.lsplit;
    cprims[numChildren] = splitter.rprims; cinfo[numChildren] = splitter.rinfo; csplit[numChildren] = splitter.rsplit;
    numChildren++;

    /*! find best child to split */
    if (numChildren < 4) 
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
                                                    &parent->alloc,parent->geom,
                                                    cprims[bestChild],cinfo[bestChild],csplit[bestChild],
                                                    _loop,this);
        return;
      }
    }

    /*! create an inner node */
    BVH4MB::Node* node = (BVH4MB::Node*) parent->bvh->alloc.malloc(threadIndex,sizeof(BVH4MB::Node),1 << BVH4MB::alignment); node->clear();
    dst = BVH4MB::Base::encodeNode(node);
    for (size_t i=0; i<numChildren; i++) {
      node->set(i,cinfo[i].geomBounds,NULL);
      parent->recurse(threadIndex,threadCount,event,node->child[i],depth+1,cprims[i],cinfo[i],csplit[i]);
    }
    delete this;
  }

  /* explicit template instantiations */
  INSTANTIATE_TEMPLATE_BY_HEURISTIC(BVH4MBBuilder);
}
