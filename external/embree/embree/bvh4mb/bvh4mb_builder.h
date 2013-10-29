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

#ifndef __EMBREE_BVH4MB_BUILDER_H__
#define __EMBREE_BVH4MB_BUILDER_H__

#include "bvh4mb.h"
#include "../builders/primrefalloc.h"
#include "../builders/primrefblock.h"
#include "../builders/primrefgen.h"
#include "../builders/splitter.h"
#include "../builders/splitter_parallel.h"

namespace embree
{
  /* BVH4MB builder. The builder is multi-threaded and implements 3
   * different build strategies: 1) Small tasks are finished in a
   * single thread (BuildTask) 2) Medium sized tasks are split into
   * two tasks using a single thread (SplitTask) and 3) Large tasks are
   * split using multiple threads on one processor. */

  template<typename Heuristic>
    class BVH4MBBuilder : public RefCount
  {
    ALIGNED_CLASS;
  public:

    /*! Type of BVH build */
    typedef BVH4MB Type;

    /*! Split type of the split heuristic. */
    typedef typename Heuristic::Split Split;
    typedef typename Heuristic::PrimInfo PrimInfo;
    
    /*! Normal splitters based on heuristic */
    typedef PrimRefGen<Heuristic,atomic_set<PrimRefBlock> > PrimRefGenNormal;
    typedef MultiThreadedSplitter<Heuristic,atomic_set<PrimRefBlock> > MultiThreadedSplitterNormal;
    typedef Splitter<Heuristic> SplitterNormal;
    
  public:

    /*! creates the acceleration structure */
    static void create (TaskScheduler::Event* event, Accel* accel, const size_t minLeafSize = 1, const size_t maxLeafSize = inf) {
      new BVH4MBBuilder(event,(BVH4MB*)accel,minLeafSize,maxLeafSize);
    }

    /*! Constructor. */
    BVH4MBBuilder(TaskScheduler::Event* event, BVH4MB* bvh, const size_t minLeafSize = 1, const size_t maxLeafSize = inf);

    /*! Destructor. */
    ~BVH4MBBuilder();

    /*! build job */
    TASK_COMPLETE_FUNCTION(BVH4MBBuilder,build);

    /*! finishes the build */
    TASK_COMPLETE_FUNCTION(BVH4MBBuilder,finish);

    /*! creates a leaf node */
    BVH4MB::Base* createLeaf(size_t threadIndex, atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo);

    /*! Selects between full build and single-threaded split strategy. */
    void recurse(size_t threadIndex, size_t threadCount, TaskScheduler::Event* group, 
                 BVH4MB::Base*& node, size_t depth, atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split);

    /*! Single-threaded task that builds a complete BVH. */
    class BuildTask {
      ALIGNED_CLASS
    public:

      /*! Default task construction. */
      BuildTask(size_t threadIndex, size_t threadCount, TaskScheduler::Event* group, 
                BVH4MBBuilder* parent, BVH4MB::Base*& node, size_t depth, 
                atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split);

      /*! Task entry function. */
      TASK_COMPLETE_FUNCTION(BuildTask,run);

      /*! Recursively finishes the BVH construction. */
      BVH4MB::Base* recurse(size_t depth, atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split);

    private:
      size_t threadIndex;
      TaskScheduler::Task task;
      BVH4MBBuilder*                   parent;  //!< Pointer to parent task.
      BVH4MB::Base*&                   dst;     //!< Reference to output the node.
      size_t                           depth;   //!< Recursion depth of the root of this subtree.
      atomic_set<PrimRefBlock>         prims;   //!< The list of primitives.
      PrimInfo                         pinfo;   //!< Bounding info of primitives.
      Split                            split;   //!< The best split for the primitives.
    };

    /*! Single-threaded task that builds a single node and creates subtasks for the children. */
    class SplitTask {
      ALIGNED_CLASS
    public:

      /*! Default task construction. */
      SplitTask(size_t threadIndex, size_t threadCount, TaskScheduler::Event* group, 
                BVH4MBBuilder* parent, BVH4MB::Base*& node, size_t depth, 
                atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split);

      /*! Task entry function. */
      TASK_COMPLETE_FUNCTION(SplitTask,recurse);

    private:
      TaskScheduler::Task task;
      BVH4MBBuilder*  parent;         //!< Pointer to parent task.
      BVH4MB::Base*&  dst;            //!< Reference to output the node.
      size_t          depth;          //!< Recursion depth of this node.
      atomic_set<PrimRefBlock> prims; //!< The list of primitives.
      PrimInfo        pinfo;          //!< Bounding info of primitives.
      Split           split;          //!< The best split for the primitives.
    };

    /*! Multi-threaded split task that builds a single node and creates subtasks for the children. */
    class ParallelSplitTask {
      ALIGNED_CLASS
    public:

      /*! Default task construction. */
      ParallelSplitTask(size_t threadIndex, size_t threadCount, TaskScheduler::Event* group, 
                        BVH4MBBuilder* parent, BVH4MB::Base*& node, size_t depth, 
                        atomic_set<PrimRefBlock>& prims, const PrimInfo& pinfo, const Split& split);

      /*! Called after the parallel binning to creates the node. */
      TASK_COMPLETE_FUNCTION(ParallelSplitTask,loop);

    private:
      BVH4MBBuilder*              parent;                 //!< Pointer to parent task.
      BVH4MB::Base*&              dst;                    //!< Reference to output the node ID.
      size_t                      depth;                  //!< Recursion depth of this node.
      atomic_set<PrimRefBlock>    cprims[4];              //!< Primitive lists for the children.
      PrimInfo                    cinfo [4];              //!< Bounding info for the children.
      Split                       csplit[4];              //!< Best next split for the children.
      size_t                      numChildren;            //!< Current number of children.
      ssize_t                     bestChild;              //!< Best child selected for splitting.
      MultiThreadedSplitterNormal splitter;               //!< Parallel splitter
    };

  private:
    RTCGeometry* geom;                //!< input geometry
    const TriangleType& trity;          //!< triangle type stored in BVH
    size_t minLeafSize;                 //!< minimal size of a leaf
    size_t maxLeafSize;                 //!< maximal size of a leaf
    PrimRefAlloc alloc;                 //!< Allocator for primitive blocks
    TaskScheduler::QUEUE taskQueue;     //!< Task queue to use

    TaskScheduler::EventScheduleTask finishStage;
    TaskScheduler::EventScheduleTask buildStage;
    PrimRefGenNormal initStage;               //!< job to generate build primitives

  public:
    BVH4MB* bvh;                    //!< Output BVH
  };
}

#endif
