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

#include "bvh2_builder.h"
#include "../common/heuristics.h"
#include "../triangle/triangles.h"
#include "sys/sync/barrier.h"
#include "sys/thread.h"

namespace embree
{
  template<typename Heuristic>
  BVH2Builder<Heuristic>::BVH2Builder(const TriangleType& trity, const std::string& intTy,
                                      const BuildTriangle* triangles, size_t numTriangles, 
                                      const Vec3fa* vertices, size_t numVertices, const BBox3f& bounds, bool freeArrays)
    : trity(trity), triangles(triangles), numTriangles(numTriangles), vertices(vertices), numVertices(numVertices), 
      taskQueue(Heuristic::depthFirst ? TaskScheduler::GLOBAL_FRONT : TaskScheduler::GLOBAL_BACK),
      bvh(new BVH2(trity,intTy,vertices,numVertices,freeArrays))
  {
    /* generate build primitives */
    scheduler->start();
    PrimRefGenNormal gen(TaskScheduler::ThreadInfo(),triangles,numTriangles,vertices,numVertices,bounds,&alloc);
    scheduler->stop();
    
    /* start parallel build */
    scheduler->start();
    recurse(TaskScheduler::ThreadInfo(),bvh->root,1,gen.prims,gen.pinfo,gen.split);
    scheduler->stop();

    /* rotate top part of tree */
    for (int i=0; i<5; i++) bvh->rotate(bvh->root,1);
    bvh->sort(bvh->root,inf);
    bvh->clearBarrier(bvh->root);
  }

  template<typename Heuristic>
  BVH2Builder<Heuristic>::~BVH2Builder() {
  }

  template<typename Heuristic>
  void BVH2Builder<Heuristic>::recurse(const TaskScheduler::ThreadInfo& thread, BVH2::Base*& node, size_t depth, 
                                       atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
  {
    /* use full single threaded build for small jobs */
    if (pinfo.size() < 4*1024)
      new BuildTask(thread,this,node,depth,prims,pinfo,split);

    /* use single threaded split for medium size jobs  */
    else if (pinfo.size() < 128*1024)
      new SplitTask(thread,this,node,depth,prims,pinfo,split);

    /* use parallel splitter for big jobs */
    else new ParallelSplitTask(thread,this,node,depth,prims,pinfo,split);
  }

  /***********************************************************************************************************************
   *                                         Leaf creation
   **********************************************************************************************************************/

  template<typename Heuristic>
  typename BVH2::Base* BVH2Builder<Heuristic>::createLeaf(const TaskScheduler::ThreadInfo& thread, atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo)
  {
    /* allocate leaf node */
    size_t blocks = trity.blocks(pinfo.size());
    char* leaf = (char*) bvh->alloc.malloc(thread,blocks*trity.bytes,1 << BVH2::alignment);

    /* insert all triangles */
    atomic_set<PrimRefBlock>::block_iterator_unsafe iter(prims);
    for (size_t i=0; i<blocks; i++) trity.pack(leaf+i*trity.bytes,iter,triangles,vertices);
    assert(!iter);

    /* free all primitive blocks */
    while (atomic_set<PrimRefBlock>::item* block = prims.take())
      alloc.free(thread,block);

    return BVH2::Base::encodeLeaf(leaf,blocks);
  }

  /***********************************************************************************************************************
   *                                         Full Recursive Build Task
   **********************************************************************************************************************/

  template<typename Heuristic>
  __forceinline BVH2Builder<Heuristic>::BuildTask::BuildTask(const TaskScheduler::ThreadInfo& thread, BVH2Builder* parent, BVH2::Base*& node, size_t depth, 
                       atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
    : thread(NULL), parent(parent), dst(node), depth(depth), prims(prims), pinfo(pinfo), split(split)
  {
    scheduler->addTask(thread,parent->taskQueue,(TaskScheduler::runFunction)&BuildTask::run,this,1,NULL,NULL,"build::full");
  }

  template<typename Heuristic>
  void BVH2Builder<Heuristic>::BuildTask::run(const TaskScheduler::ThreadInfo& thread, BuildTask* This, size_t elts)
  {
    This->thread = &thread;
    This->depth = max(This->depth,BVH2::maxDepth-BVH2::maxLocalDepth+1);
    This->dst = This->recurse(This->depth,This->prims,This->pinfo,This->split);
    for (int i=0; i<5; i++) This->parent->bvh->rotate(This->dst,This->depth);
    This->parent->bvh->sort(This->dst,inf);
    This->dst = This->dst->setBarrier();
    delete This;
  }

  template<typename Heuristic>
  BVH2::Base* BVH2Builder<Heuristic>::BuildTask::recurse(size_t depth, atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
  {
    /*! compute leaf and split cost */
    const float leafSAH  = parent->trity.intCost*pinfo.sah();
    const float splitSAH = BVH2::travCost*halfArea(pinfo.geomBounds)+parent->trity.intCost*split.sah();
    assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(prims).size() == pinfo.size());
    assert(pinfo.size() == 0 || leafSAH >= 0 && splitSAH >= 0);

    /*! create a leaf node when threshold reached or SAH tells us to stop */
    if (pinfo.size() <= 1 || depth > BVH2::maxDepth || (pinfo.size() <= parent->bvh->maxLeafTris && leafSAH <= splitSAH)) {
      return parent->createLeaf(*thread,prims,pinfo);
    }

    /*! perform best found split and find new splits */
    SplitterNormal splitter(*thread,&parent->alloc,parent->triangles,parent->vertices,prims,pinfo,split);

    /*! create an inner node */
    BVH2::Node* node = (BVH2::Node*) parent->bvh->alloc.malloc(*thread,sizeof(BVH2::Node),1 << BVH2::alignment); node->clear();
    node->set(1,splitter.rinfo.geomBounds,recurse(depth+1,splitter.rprims,splitter.rinfo,splitter.rsplit));
    node->set(0,splitter.linfo.geomBounds,recurse(depth+1,splitter.lprims,splitter.linfo,splitter.lsplit));
    return BVH2::Base::encodeNode(node);
  }

  /***********************************************************************************************************************
   *                                      Single Threaded Split Task
   **********************************************************************************************************************/

  template<typename Heuristic>
  __forceinline BVH2Builder<Heuristic>::SplitTask::SplitTask(const TaskScheduler::ThreadInfo& thread, BVH2Builder* parent, BVH2::Base*& node, size_t depth, 
                       atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
    : parent(parent), dst(node), depth(depth), prims(prims), pinfo(pinfo), split(split)
  {
    scheduler->addTask(thread,parent->taskQueue,(TaskScheduler::runFunction)_recurse,this,1,NULL,NULL,"build::split");
  }

  template<typename Heuristic>
  void BVH2Builder<Heuristic>::SplitTask::recurse(const TaskScheduler::ThreadInfo& thread)
  {
     /*! compute leaf and split cost */
    const float leafSAH  = parent->trity.intCost*pinfo.sah();
    const float splitSAH = BVH2::travCost*halfArea(pinfo.geomBounds)+parent->trity.intCost*split.sah();
    assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(prims).size() == pinfo.size());
    assert(pinfo.size() == 0 || leafSAH >= 0 && splitSAH >= 0);

    /*! create a leaf node when threshold reached or SAH tells us to stop */
    if (pinfo.size() <= 1 || depth > BVH2::maxDepth || (pinfo.size() <= parent->bvh->maxLeafTris && leafSAH <= splitSAH)) {
      dst = parent->createLeaf(thread,prims,pinfo); delete this; return;
    }

    /*! perform best found split and find new splits */
    SplitterNormal splitter(thread,&parent->alloc,parent->triangles,parent->vertices,prims,pinfo,split);

    /*! create an inner node */
    BVH2::Node* node = (BVH2::Node*) parent->bvh->alloc.malloc(thread,sizeof(BVH2::Node),1 << BVH2::alignment); node->clear();
    dst = BVH2::Base::encodeNode(node);
    node->set(1,splitter.rinfo.geomBounds,0);
    node->set(0,splitter.linfo.geomBounds,0);
    parent->recurse(thread,node->child[0],depth+1,splitter.lprims,splitter.linfo,splitter.lsplit);
    parent->recurse(thread,node->child[1],depth+1,splitter.rprims,splitter.rinfo,splitter.rsplit);
    delete this;
  }
  
  /***********************************************************************************************************************
   *                                           Parallel Split Task
   **********************************************************************************************************************/

  template<typename Heuristic>
  __forceinline BVH2Builder<Heuristic>::ParallelSplitTask::ParallelSplitTask(const TaskScheduler::ThreadInfo& thread, BVH2Builder* parent, BVH2::Base*& node, size_t depth, 
                                                                             atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
    : parent(parent), dst(node), depth(depth), split(split)
  {
    /*! compute leaf and split cost */
    const float leafSAH  = parent->trity.intCost*pinfo.sah();
    const float splitSAH = BVH2::travCost*halfArea(pinfo.geomBounds)+parent->trity.intCost*split.sah();
    assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(prims).size() == pinfo.size());
    assert(pinfo.size() == 0 || leafSAH >= 0 && splitSAH >= 0);

    /*! create a leaf node when threshold reached or SAH tells us to stop */
    if (pinfo.size() <= 1 || depth > BVH2::maxDepth || (pinfo.size() <= parent->bvh->maxLeafTris && leafSAH <= splitSAH)) {
      dst = parent->createLeaf(thread,prims,pinfo); delete this; return;
    }

    /*! start parallel splitting */
    new (&splitter) MultiThreadedSplitterNormal(thread,
                                                &parent->alloc,parent->triangles,parent->vertices,prims,pinfo,split,
                                                (TaskScheduler::completeFunction)_createNode,this);
  }

  template<typename Heuristic>
  void BVH2Builder<Heuristic>::ParallelSplitTask::createNode(const TaskScheduler::ThreadInfo& thread)
  {
    /*! create an inner node */
    BVH2::Node* node = (BVH2::Node*) parent->bvh->alloc.malloc(thread,sizeof(BVH2::Node),1 << BVH2::alignment); node->clear();
    dst = BVH2::Base::encodeNode(node);
    node->set(1,splitter.rinfo.geomBounds,0);
    node->set(0,splitter.linfo.geomBounds,0);
    parent->recurse(thread,node->child[0],depth+1,splitter.lprims,splitter.linfo,splitter.lsplit);
    parent->recurse(thread,node->child[1],depth+1,splitter.rprims,splitter.rinfo,splitter.rsplit);
    delete this;
  }
  
  /* explicit template instantiations */
  INSTANTIATE_TEMPLATE_BY_HEURISTIC(BVH2Builder);
}

