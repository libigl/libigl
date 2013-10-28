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

#ifndef __EMBREE_BVH_BUILDER_TEMPLATED_H__
#define __EMBREE_BVH_BUILDER_TEMPLATED_H__

#include "primrefalloc.h"
#include "primrefblock.h"
#include "primrefgen.h"
#include "splitter.h"
#include "splitter_fallback.h"
#include "splitter_parallel.h"
#include "../geometry/triangle.h"
#include "bvh_sort.h"
#include "bvh_rotate.h"

#define ROTATE_TREE 1
#define SORT_TREE 1

namespace embree
{
  /* BVH builder. The builder is multi-threaded and implements 3
   * different build strategies: 1) Small tasks are finished in a
   * single thread (BuildTask) 2) Medium sized tasks are split into
   * two tasks using a single thread (SplitTask) and 3) Large tasks are
   * split using multiple threads on one processor. */

  template<typename BVH, typename Heuristic>
    class BVHBuilderT : public RefCount
  {
    ALIGNED_CLASS;
  public:

    /*! Type shortcuts */
    typedef BVH Type;
    typedef typename BVH::Node    Node;
    typedef typename BVH::NodeRef NodeRef;

    /*! Split type of the split heuristic. */
    typedef typename Heuristic::Split Split;
    typedef typename Heuristic::PrimInfo PrimInfo;
    
    /*! Normal splitters based on heuristic */
    typedef PrimRefGen<Heuristic,atomic_set<PrimRefBlock> > PrimRefGenNormal;
    typedef MultiThreadedSplitter<Heuristic,atomic_set<PrimRefBlock> > MultiThreadedSplitterNormal;
    typedef Splitter<Heuristic> SplitterNormal;
    
  public:

    /*! creates the acceleration structure */
    static void create (TaskScheduler::Event* event, Accel* accel, const size_t minLeafSize = inf, const size_t maxLeafSize = inf) { 
      new BVHBuilderT(event,(BVH*)accel,minLeafSize,maxLeafSize);
    }

    /*! Constructor. */
    BVHBuilderT (TaskScheduler::Event* event, BVH* bvh, const size_t minLeafSize = 1, const size_t maxLeafSize = inf)
      : geom(bvh->geom), trity(bvh->trity),
        minLeafSize(minLeafSize), maxLeafSize(maxLeafSize),
        taskQueue(Heuristic::depthFirst ? TaskScheduler::GLOBAL_BACK : TaskScheduler::GLOBAL_FRONT),
        bvh(bvh)
    {
      if (bvh->maxLeafTris < this->maxLeafSize) 
        this->maxLeafSize = bvh->maxLeafTris;

      /* calculate required memory */
      size_t tris  = bvh->geom->size();
      size_t prims = 2*bvh->trity.blocks(tris); // reserve more memory than typically needed
      size_t nodes = 4*(prims+BVH::N-1)/BVH::N; // reserve more memory than typically needed
      bvh->init(nodes,prims);

      /* create tasks */
      new (&finishStage) TaskScheduler::EventScheduleTask(event,_finish,this,"BVHBuilder::finish");
      new (&buildStage ) TaskScheduler::EventScheduleTask(&finishStage,_build,this,"BVHBuilder::build");
      new (&initStage  ) PrimRefGenNormal(&buildStage,geom,&alloc);
    }

    /*! build job */
    TASK_COMPLETE_FUNCTION_(BVHBuilderT,build);
    void build(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event) {
      recurse(threadIndex,threadCount,event,bvh->root,1,initStage.prims,initStage.pinfo,initStage.split);
    }

    /*! finishes the build */
    TASK_COMPLETE_FUNCTION_(BVHBuilderT,finish);
    void finish(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event) 
    {
#if ROTATE_TREE
      for (int i=0; i<5; i++) BVHRotateT<BVH>::rotate(bvh,bvh->root,1,inf);
#endif
#if SORT_TREE
      BVHSortT<BVH>::sort(geom,bvh,bvh->root);
#endif
      bvh->clearBarrier(bvh->root);
      bvh->bounds = initStage.pinfo.geomBounds;
      delete this;
    }

    /*! creates a leaf node */
    NodeRef createLeafNode(size_t threadIndex, atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo)
    {
      /* allocate leaf node */
      size_t blocks = trity.blocks(pinfo.size());
      char* leaf = bvh->allocTris(threadIndex,blocks);
      
      /* insert all triangles */
      atomic_set<PrimRefBlock>::block_iterator_unsafe iter(prims);
      for (size_t i=0; i<blocks; i++) {
        trity.pack(leaf+i*trity.bytes,iter,geom);
      }
      assert(!iter);
      
      /* free all primitive blocks */
      while (atomic_set<PrimRefBlock>::item* block = prims.take())
        alloc.free(threadIndex,block);
      
      return BVH::NodeRef::encodeLeaf(bvh->triPtr(),leaf,blocks);
    }

    NodeRef createLeaf(size_t threadIndex, atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, size_t depth)
    {
#if defined(_DEBUG)
      if (depth >= BVH::maxBuildDepthLeaves) 
        throw std::runtime_error("ERROR: BVH too deep.");
#endif
      
      /* create leaf for few primitives */
      if (pinfo.size() <= maxLeafSize)
        return createLeafNode(threadIndex,prims,pinfo);
      
      /* perform split */
      atomic_set<PrimRefBlock> cprims[2];
      PrimInfo                 cinfo[2];
      FallBackSplitter<Heuristic,atomic_set<PrimRefBlock> >::split(threadIndex,&alloc,geom,prims,pinfo,cprims[0],cinfo[0],cprims[1],cinfo[1]);
      
      /*! create an inner node */
      Node* node = bvh->allocNode(threadIndex);
      for (size_t i=0; i<2; i++) node->set(i,cinfo[i].geomBounds,createLeaf(threadIndex,cprims[i],cinfo[i],depth+1));
      return BVH::NodeRef::encodeNode(bvh->nodePtr(),node);
    }  

    /*! Selects between full build and single-threaded split strategy. */
    void recurse(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event, 
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

    /***********************************************************************************************************************
     *                                      Single Threaded Build Task
     **********************************************************************************************************************/
    
    /*! Single-threaded task that builds a complete BVH. */
    class BuildTask {
      ALIGNED_CLASS
    public:

      /*! Default task construction. */
      BuildTask(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event, 
                BVHBuilderT* parent, NodeRef& node, size_t depth, 
                atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
        : threadIndex(threadIndex), parent(parent), dst(node), depth(depth), prims(prims), pinfo(pinfo), split(split)
      {
        new (&task) TaskScheduler::Task(event,_run,this,"build::full");
        TaskScheduler::addTask(threadIndex,parent->taskQueue,&task);
      }

      /*! Task entry function. */
      TASK_COMPLETE_FUNCTION_(BuildTask,run);
      void run(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
      {
        this->threadIndex = threadIndex;
        dst = recurse(depth,prims,pinfo,split);
#if ROTATE_TREE
        for (int i=0; i<5; i++) BVHRotateT<BVH>::rotate(parent->bvh,dst,depth,inf); 
#endif
#if SORT_TREE
        BVHSortT<BVH>::sort(parent->geom,parent->bvh,dst); 
#endif
        dst.setBarrier();
        delete this;
      }

      /*! Recursively finishes the BVH construction. */
      NodeRef recurse(size_t depth, atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
      {
        /*! compute leaf and split cost */
        const float leafSAH  = parent->trity.intCost*pinfo.sah();
        const float splitSAH = BVH::travCost*halfArea(pinfo.geomBounds)+parent->trity.intCost*split.sah();
        assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(prims).size() == pinfo.size());
        assert(pinfo.size() == 0 || leafSAH >= 0 && splitSAH >= 0);
        
        /*! create a leaf node when threshold reached or SAH tells us to stop */
        if (pinfo.size() <= parent->minLeafSize || depth > BVH::maxBuildDepth || (pinfo.size() <= parent->maxLeafSize && leafSAH <= splitSAH)) {
          return parent->createLeaf(threadIndex,prims,pinfo,depth);
        }
        
        /*! initialize child list */
        atomic_set<PrimRefBlock> cprims[BVH::N]; cprims[0] = prims;
        PrimInfo                 cinfo [BVH::N]; cinfo [0] = pinfo;
        Split                    csplit[BVH::N]; csplit[0] = split;
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
          
        } while (numChildren < BVH::N);
        
        /*! create an inner node */
        Node* node = parent->bvh->allocNode(threadIndex);
        for (size_t i=0; i<numChildren; i++) node->set(i,cinfo[i].geomBounds,recurse(depth+1,cprims[i],cinfo[i],csplit[i]));
        return BVH::NodeRef::encodeNode(parent->bvh->nodePtr(),node);
      }
      
    private:
      size_t threadIndex;
      TaskScheduler::Task task;
      
      BVHBuilderT*                     parent;   //!< Pointer to parent task.
      NodeRef&                     dst;      //!< Reference to output the node.
      size_t                           depth;    //!< Recursion depth of the root of this subtree.
      atomic_set<PrimRefBlock>         prims;    //!< The list of primitives.
      PrimInfo                         pinfo;    //!< Bounding info of primitives.
      Split                            split;    //!< The best split for the primitives.
    };

    /***********************************************************************************************************************
     *                                      Single Threaded Split Task
     **********************************************************************************************************************/

    /*! Single-threaded task that builds a single node and creates subtasks for the children. */
    class SplitTask {
      ALIGNED_CLASS
    public:

      /*! Default task construction. */
      SplitTask(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event, 
                BVHBuilderT* parent, NodeRef& node, size_t depth, 
                atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
        : parent(parent), dst(node), depth(depth), prims(prims), pinfo(pinfo), split(split)
      {
        new (&task) TaskScheduler::Task(event,_recurse,this,"build::split");
        TaskScheduler::addTask(threadIndex,parent->taskQueue,&task);
      }

      /*! Task entry function. */
      TASK_COMPLETE_FUNCTION_(SplitTask,recurse);
      void recurse(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
      {
        /*! compute leaf and split cost */
        const float leafSAH  = parent->trity.intCost*pinfo.sah();
        const float splitSAH = BVH::travCost*halfArea(pinfo.geomBounds)+parent->trity.intCost*split.sah();
        assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(prims).size() == pinfo.size());
        assert(pinfo.size() == 0 || leafSAH >= 0 && splitSAH >= 0);
        
        /*! create a leaf node when threshold reached or SAH tells us to stop */
        if (pinfo.size() <= parent->minLeafSize || depth > BVH::maxBuildDepth || (pinfo.size() <= parent->maxLeafSize && leafSAH <= splitSAH)) {
          dst = parent->createLeaf(threadIndex,prims,pinfo,depth); delete this; return;
        }
        
        /*! initialize child list */
        atomic_set<PrimRefBlock> cprims[BVH::N]; 
        PrimInfo                 cinfo [BVH::N]; 
        Split                    csplit[BVH::N]; 
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
          SplitterNormal splitter(threadIndex,&parent->alloc,parent->geom,cprims[bestChild],cinfo[bestChild],csplit[bestChild]);
          cprims[bestChild  ] = splitter.lprims; cinfo[bestChild  ] = splitter.linfo; csplit[bestChild  ] = splitter.lsplit;
          cprims[numChildren] = splitter.rprims; cinfo[numChildren] = splitter.rinfo; csplit[numChildren] = splitter.rsplit;
          numChildren++;
          
        } while (numChildren < BVH::N);
        
        /*! create an inner node */
        Node* node = parent->bvh->allocNode(threadIndex);
        dst = BVH::NodeRef::encodeNode(parent->bvh->nodePtr(),node);
        for (size_t i=0; i<numChildren; i++) {
          node->set(i,cinfo[i].geomBounds,0);
          parent->recurse(threadIndex,threadCount,event,node->child(i),depth+1,cprims[i],cinfo[i],csplit[i]);
        }
        delete this;
      }
      
    private:
      TaskScheduler::Task task;
      BVHBuilderT*    parent;         //!< Pointer to parent task.
      NodeRef&    dst;            //!< Reference to output the node.
      size_t          depth;          //!< Recursion depth of this node.
      atomic_set<PrimRefBlock> prims; //!< The list of primitives.
      PrimInfo        pinfo;          //!< Bounding info of primitives.
      Split           split;          //!< The best split for the primitives.
    };

  /***********************************************************************************************************************
   *                                           Parallel Split Task
   **********************************************************************************************************************/

    /*! Multi-threaded split task that builds a single node and creates subtasks for the children. */
    class ParallelSplitTask {
      ALIGNED_CLASS
    public:

      /*! Default task construction. */
      ParallelSplitTask (size_t threadIndex, size_t threadCount, TaskScheduler::Event* event, 
                         BVHBuilderT* parent, NodeRef& node, size_t depth, 
                         atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split)
        : parent(parent), dst(node), depth(depth), numChildren(1), bestChild(0)
      {
        /*! compute leaf and split cost */
        const float leafSAH  = parent->trity.intCost*pinfo.sah();
        const float splitSAH = BVH::travCost*halfArea(pinfo.geomBounds)+parent->trity.intCost*split.sah();
        assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(prims).size() == pinfo.size());
        assert(pinfo.size() == 0 || leafSAH >= 0 && splitSAH >= 0);
        
        /*! create a leaf node when threshold reached or SAH tells us to stop */
        if (pinfo.size() <= parent->minLeafSize || depth > BVH::maxBuildDepth || (pinfo.size() <= parent->maxLeafSize && leafSAH <= splitSAH)) {
          dst = parent->createLeaf(threadIndex,prims,pinfo,depth); delete this; return;
        }
        
        /*! initialize child list */
        cprims[0] = prims; cinfo [0] = pinfo; csplit[0] = split; 
        
        /*! perform first split */
        new (&splitter) MultiThreadedSplitterNormal(threadIndex,threadCount,event,
                                                    &parent->alloc,parent->geom,
                                                    cprims[bestChild],cinfo[bestChild],csplit[bestChild],
                                                    _loop,this);
      }

      /*! Called after the parallel binning to creates the node. */
      TASK_COMPLETE_FUNCTION_(ParallelSplitTask,loop);
      void loop(size_t threadIndex, size_t threadCount, TaskScheduler::Event* event)
      {
        /*! copy new children into work array */
        cprims[bestChild  ] = splitter.lprims; cinfo[bestChild  ] = splitter.linfo; csplit[bestChild  ] = splitter.lsplit;
        cprims[numChildren] = splitter.rprims; cinfo[numChildren] = splitter.rinfo; csplit[numChildren] = splitter.rsplit;
        numChildren++;
        
        /*! find best child to split */
        if (numChildren < BVH::N) 
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
        Node* node = parent->bvh->allocNode(threadIndex);
        dst = BVH::NodeRef::encodeNode(parent->bvh->nodePtr(),node);
        for (size_t i=0; i<numChildren; i++) {
          node->set(i,cinfo[i].geomBounds,0);
          parent->recurse(threadIndex,threadCount,event,node->child(i),depth+1,cprims[i],cinfo[i],csplit[i]);
        }
        delete this;
      }
      
    private:
      BVHBuilderT*                parent;                 //!< Pointer to parent task.
      NodeRef&                    dst;                    //!< Reference to output the node ID.
      size_t                      depth;                  //!< Recursion depth of this node.
      atomic_set<PrimRefBlock>    cprims[BVH::N];              //!< Primitive lists for the children.
      PrimInfo                    cinfo [BVH::N];              //!< Bounding info for the children.
      Split                       csplit[BVH::N];              //!< Best next split for the children.
      size_t                      numChildren;            //!< Current number of children.
      ssize_t                     bestChild;              //!< Best child selected for splitting.
      MultiThreadedSplitterNormal splitter;               //!< Parallel splitter
    };

  private:
    RTCGeometry* geom;      //!< input geometry

  public:
    const TriangleType& trity;          //!< triangle type stored in BVH
    size_t minLeafSize;                 //!< minimal size of a leaf
    size_t maxLeafSize;                 //!< maximal size of a leaf
    PrimRefAlloc alloc;                 //!< Allocator for primitive blocks
    TaskScheduler::QUEUE taskQueue;     //!< Task queue to use

    TaskScheduler::EventScheduleTask finishStage;
    TaskScheduler::EventScheduleTask buildStage;
    PrimRefGenNormal initStage;               //!< job to generate build primitives

  public:
    BVH* bvh;                      //!< Output BVH
  };
}

#endif
