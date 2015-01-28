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

#include "bvh4i.h"
#include "bvh4i_builder_morton_enhanced.h"
#include "bvh4i_builder_binner.h"
#include "bvh4i_statistics.h"

#include "bvh4i_builder_util.h"
#include "limits.h"
#include "sys/sync/barrier.h"

#define TREE_BRANCHING_FACTOR 4
#define QBVH_BUILDER_LEAF_ITEM_THRESHOLD 4

#define DIAG_FACTOR 0.2f
#define MAX_REBUILD_NODES 1024*2

//#define PROFILE
//#define TEST_BUILD_PERFORMANCE

#define DBG(x) 

namespace embree 
{
  namespace isa
  {
    BVH4iBuilderMortonEnhanced::BVH4iBuilderMortonEnhanced (BVH4i* _bvh, BuildSource* _source, void* _geometry, const size_t _minLeafSize, const size_t _maxLeafSize)
    : BVH4iBuilderMorton(_bvh,_source,_geometry,_minLeafSize,_maxLeafSize){}
    
    void BVH4iBuilderMortonEnhanced::build(size_t threadIndex, size_t threadCount) 
    {
      init();
      TaskScheduler::executeTask(threadIndex,threadCount,_build_parallel_morton_enhanced,this,TaskScheduler::getNumThreads(),"build_parallel_morton_enhanced");
    }
    
    bool splitSAH(PrimRef * __restrict__ const primref, BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild)
    {
      /* mark as leaf if leaf threshold reached */
      const unsigned int items = current.end - current.begin;
      if (items <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD) {
        current.createLeaf();
        return false;
      }
      
      /* calculate binning function */
      Mapping<16> mapping(current.bounds);
      
      /* binning of centroids */
      Binner<16> binner;
      binner.bin(primref,current.begin,current.end,mapping);
      
      /* find best split */
      Split split; 
      binner.best(split,mapping);
      
      /* if we cannot find a valid split, enforce an arbitrary split */
      if (unlikely(split.pos == -1)) 
        return split_fallback(primref,current,leftChild,rightChild);
      
      /* partitioning of items */
      binner.partition(primref, current.begin, current.end, split, mapping, leftChild, rightChild);
      if (leftChild.items()  <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD) leftChild.createLeaf();
      if (rightChild.items() <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD) rightChild.createLeaf();	
      return true;
    }
    
    
    
    void BVH4iBuilderMortonEnhanced::extractTopLevelTree(const size_t index,
                                                         const Vec3fa &root_diag,
                                                         BVHNode *__restrict__ const local_node,
                                                         size_t &nodes)
    {
      BVHNode &entry = this->node[index];
      
      if (!entry.isLeaf())
      {
	const unsigned int children = entry.firstChildID();
        
	const Vec3fa diag = entry.upper - entry.lower;
	if (diag.x < DIAG_FACTOR * root_diag.x &&
	    diag.y < DIAG_FACTOR * root_diag.y &&
	    diag.z < DIAG_FACTOR * root_diag.z)
        {
          if (unlikely(nodes >= MAX_REBUILD_NODES)) FATAL("too many subtrees");
          //subTreeIDs[subTrees++] = index;	      	      
          local_node[nodes++] = entry;
          assert(nodes < MAX_REBUILD_NODES);
          return;
        }
	else
        {
          for (unsigned int i=0;i<entry.items();i++) 
          {
            const unsigned int childIndex = children + i;	    
            extractTopLevelTree(childIndex,root_diag,local_node,nodes);	  
          }
        }
      }
      else
      {
	if (unlikely(nodes >= MAX_REBUILD_NODES)) FATAL("too many subtrees");
	//subTreeIDs[subTrees++] = index;	      	      
	local_node[nodes++] = entry;
	assert(nodes < MAX_REBUILD_NODES);
      }     
    }
    
    
    void BVH4iBuilderMortonEnhanced::buildTopLevelSAHTree(BVHNode &parent,
                                                          BuildRecord &current,
                                                          BVHNode *__restrict__ local_node)
    {
      DBG(std::cout << std::endl; PING);
      DBG(DBG_PRINT(parent));
      DBG(DBG_PRINT(current));
#ifdef DEBUG
      {
        PrimRef tmp;
        tmp.lower = current.bounds.geometry.lower;//storeSceneAABB((float*)&tmp);
        tmp.upper = current.bounds.geometry.upper;
        for (size_t i=0;i<current.items();i++)
          assert(subset(*(PrimRef*)&local_node[current.begin+i],tmp) == true);
        // for (size_t i=0;i<current.items();i++)
        //   cout << i << " " << local_node[current.begin+i] << endl;
      }
#endif      
      
      BuildRecord record[TREE_BRANCHING_FACTOR];
      BuildRecord left, right;
      
      unsigned int numSplits = 1;
      record[0] = current;
      
      while(1)
      {
	bool couldSplit = false;
	if (numSplits < TREE_BRANCHING_FACTOR)
        {
          int index = -1;
          float maxArea = -1.0f;
          for (unsigned int i=0;i<numSplits;i++)
          {
            if (!record[i].isLeaf())
            {
              assert(record[i].sceneArea() > 0.0f);
              
              if (record[i].sceneArea() >= maxArea)
              {
                maxArea = record[i].sceneArea();
                index = i;
              }
            }
          }
          
          assert(index < TREE_BRANCHING_FACTOR);
          if (index != -1)
          {
            bool s = splitSAH((PrimRef*)local_node,record[index],left,right);  
            
            assert(numSplits > 0);
            
            if ( s )
            {
              record[index] = record[numSplits-1];
              numSplits--;
              
              record[numSplits+0] = left;
              record[numSplits+1] = right;
              numSplits+=2;
              
              couldSplit = true;
            }
            else 
            {// became leaf
              DBG(std::cout << "became leaf" << std::endl);
              continue;
            }
          }
        }
	// ==================================================
	if (!couldSplit) break;
      }
      
      DBG(std::cout << "DETERMINED ALL SPLITS" << std::endl);
      
      
      DBG(DBG_PRINT(numSplits));
      
      // could not generate any split, propagate leaf back to parent
      if ( numSplits == 1 )
      {
	current = record[0];
	DBG(std::cout << "CREATING LOCAL NODE LEAF WITH " << current.items() << " NODES" << std::endl);
	assert(current.isLeaf());
	if(current.items() > QBVH_BUILDER_LEAF_ITEM_THRESHOLD) 
        {
          std::cout << "WARNING UNSPLITABLE LEAF " << std::endl;
          FATAL("SHOULD NOT HAPPEN FOR TOP-LEVEL TREE");
          DBG_PRINT(current.items());
          record[1] = record[0];
          record[0].end   = record[0].begin + 4;
          record[1].begin = record[0].begin + 4;
          numSplits = 2;
        }
	else
        {
          if (unlikely(current.items() == 1))
          {
            parent = local_node[current.begin];
            return;
          }
          assert(current.items() > 1);
          
          const unsigned int currentIndex = this->atomicID.add(TREE_BRANCHING_FACTOR);
          if (unlikely(currentIndex >= this->numAllocatedNodes))
          {
            DBG_PRINT(this->numAllocatedNodes);
            FATAL("not enough nodes allocated");
          }
          
          DBG(DBG_PRINT(currentIndex));
          
          
          const avxf init_node = load8f((float*)initQBVHNode);
          
#pragma unroll(TREE_BRANCHING_FACTOR)
          for (size_t i=0;i<TREE_BRANCHING_FACTOR;i++)
            store8f((float*)&this->node[currentIndex+i],init_node);
          
          parent.createNode(currentIndex,current.items());
          //current.childrenID = currentIndex;
          
          for (size_t i=0;i<current.items();i++)
          {
            assert(subset(*(PrimRef*)&local_node[current.begin+i],*(PrimRef*)&parent));
            this->node[currentIndex+i] = local_node[current.begin+i];
          }
          
          DBG(DBG_PRINT(parent));
          return;
        }
      }
      
      assert(numSplits >= 2 && numSplits <= TREE_BRANCHING_FACTOR);
      
      // ==== aquire next four nodes ====
      const unsigned int currentIndex = this->atomicID.add(TREE_BRANCHING_FACTOR);
      if (unlikely(currentIndex >= this->numAllocatedNodes))
      {
	DBG_PRINT(this->numAllocatedNodes);
	FATAL("not enough nodes allocated");
      }
      
      parent.createNode(currentIndex,numSplits);
      //current.childrenID = currentIndex;
      //assert(current.childrenID > 0);
      
      // ==== init all four nodes ====
      const avxf init_node = load8f((float*)initQBVHNode);
      
      for (size_t i=numSplits;i<TREE_BRANCHING_FACTOR;i++)
        store8f_nt((float*)&this->node[currentIndex+i],init_node);
      
      //node[current.parentID].createNode(current.childrenID,numSplits);
      
      
      for (unsigned int i=0;i<numSplits;i++)
      {
	DBG(DBG_PRINT(currentIndex+i));
	//record[i].bounds.storeSceneAABB((float*)&this->node[currentIndex+i]);
        this->node[currentIndex+i].lower = record[i].bounds.geometry.lower;
        this->node[currentIndex+i].upper = record[i].bounds.geometry.upper;
	buildTopLevelSAHTree(this->node[currentIndex+i],record[i],local_node);
      }
      
    }
    
    
    void BVH4iBuilderMortonEnhanced::build_parallel_morton_enhanced(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event) 
    {
      initThreadState(threadIndex,threadCount);
      
      if (threadIndex == 0)
      {
	LockStepTaskScheduler::init(TaskScheduler::getNumThreads()); 
        
	init();
      }
      
      barrier.wait(threadIndex,threadCount);
      
      if (threadIndex == 0)
      {
        
        
#ifdef PROFILE
	while(1)
#endif
        {
          // ===================          
          // === start timer ===
          // ===================
          
          double t0 = 0.0f;
          if (g_verbose >= 2) {
            std::cout << "building BVH4i with Enhanced Morton builder ... " << std::endl << std::flush;
            t0 = getSeconds();
          }
          
          build_main(threadIndex,taskCount);
          
          bvh->accel = this->accel;
          bvh->qbvh  = this->node;
          
          // ==============================        
          // === rebuild top level tree ===
          // ==============================
          
          const Vec3fa rootDiag = this->node[0].upper - this->node[0].lower;
          
          __aligned(64) BVHNode local_node[MAX_REBUILD_NODES];
          size_t nodes = 0;
          extractTopLevelTree(0,rootDiag,local_node,nodes);
          BuildRecord topLevelBuildRecord;
          topLevelBuildRecord.init(this->global_bounds,0,nodes);
          buildTopLevelSAHTree(this->node[0],topLevelBuildRecord,local_node);
	  
#ifdef DEBUG
          // std::cout << "starting slow tree check..." << std::endl << std::flush;
          // checkBVH4iTree(this->node,this->accel,this->numPrimitives,true);
          // std::cout << "check done" << std::endl << std::flush;
#endif
          
          // ===================================        
          // === convert to optimized layout ===
          // ===================================
          
          this->numNodes = this->atomicID >> 2;
          
          LockStepTaskScheduler::dispatchTask( task_convertToSOALayout, this, threadIndex, threadCount );
          
          
          const QBVHNode      *const __restrict__ qbvh  = (QBVHNode*)bvh->qbvh;
          bvh->root = qbvh[0].min_d[0]; 
          bvh->bounds = BBox3f(Vec3fa(qbvh->min_x[0],qbvh->min_y[0],qbvh->min_y[0]),
                               Vec3fa(qbvh->max_x[0],qbvh->max_y[0],qbvh->max_y[0]));
          
          // ==================          
          // === stop timer ===
          // ==================
          
          if (g_verbose >= 2) {
            double t1 = getSeconds();
            double perf = source->size()/(t1-t0)*1E-6;
            std::cout << "[DONE] " << t1-t0 << "sec (" << perf << " Mtris/s)" << std::endl;
            //
          }
        }
	//freeData();
	if (g_verbose >= 2) 
	  std::cout << BVH4iStatistics(bvh).str();
        
	LockStepTaskScheduler::releaseThreads(threadCount);
      }
      else
        LockStepTaskScheduler::dispatchTaskMainLoop(threadIndex,threadCount);
    };

    Builder* BVH4iTriangle1BuilderMortonEnhanced (void* bvh, BuildSource* source, Scene* scene, const size_t minLeafSize, const size_t maxLeafSize) {
      return new BVH4iBuilderMortonEnhanced((BVH4i*)bvh,source,scene,minLeafSize,maxLeafSize);
    }
  } 
}
