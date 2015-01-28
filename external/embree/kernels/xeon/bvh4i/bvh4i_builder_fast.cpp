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
#include "bvh4i_builder_fast.h"
#include "bvh4i_statistics.h"
#include "bvh4i_builder_util.h"
#include "bvh4i_builder_binner.h"

#define BVH_NODE_PREALLOC_FACTOR 1.0f
#define QBVH_BUILDER_LEAF_ITEM_THRESHOLD 4
#define THRESHOLD_FOR_SUBTREE_RECURSION 128
#define BUILD_RECORD_SPLIT_THRESHOLD 512

#define DBG(x) 

//#define PROFILE

namespace embree
{
  namespace isa
  {
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    
    static double dt = 0.0f;
    
    BVH4iBuilderFast::BVH4iBuilderFast (BVH4i* bvh, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize)
    : source(source), geometry(geometry), primTy(bvh->primTy), bvh(bvh), numPrimitives(0), numNodes(0), numAllocatedNodes(0), prims(NULL), node(NULL), accel(NULL) {}
    
    void BVH4iBuilderFast::build(size_t threadIndex, size_t threadCount) 
    {
      if (g_verbose >= 1)
        std::cout << "building BVH4i with fast SAH builder ... " << std::flush;
      
      /* allocate BVH data */
      allocateData();
      
      LockStepTaskScheduler::init(TaskScheduler::getNumThreads());
      
#if defined(PROFILE)
      PING;
      double dt_min = pos_inf;
      double dt_avg = 0.0f;
      double dt_max = neg_inf;
      for (size_t i=0; i<20; i++) 
      {
        TaskScheduler::executeTask(threadIndex,threadCount,_build_parallel,this,TaskScheduler::getNumThreads(),"build_parallel");
        dt_min = min(dt_min,dt);
        dt_avg = dt_avg + dt;
        dt_max = max(dt_max,dt);
      }
      dt_avg /= double(20);
      
      std::cout << "[DONE]" << std::endl;
      std::cout << "  min = " << 1000.0f*dt_min << "ms (" << source->size()/dt_min*1E-6 << " Mtris/s)" << std::endl;
      std::cout << "  avg = " << 1000.0f*dt_avg << "ms (" << source->size()/dt_avg*1E-6 << " Mtris/s)" << std::endl;
      std::cout << "  max = " << 1000.0f*dt_max << "ms (" << source->size()/dt_max*1E-6 << " Mtris/s)" << std::endl;
      std::cout << BVH4iStatistics(bvh).str();
      
#else
      
      TaskScheduler::executeTask(threadIndex,threadCount,_build_parallel,this,TaskScheduler::getNumThreads(),"build_parallel");
      
      if (g_verbose >= 2) {
        double perf = source->size()/dt*1E-6;
        std::cout << "[DONE] " << 1000.0f*dt << "ms (" << perf << " Mtris/s)" << std::endl;
        std::cout << BVH4iStatistics(bvh).str();
      }
#endif
      
    }
    
    void BVH4iBuilderFast::allocateData()
    {
      size_t numPrimitivesOld = numPrimitives;
      numPrimitives = source->size();
      
      const size_t additional_size = 16 * CACHELINE_SIZE;
      const size_t numPrimitives = this->numPrimitives;
      
      if (numPrimitivesOld != numPrimitives)
      {
        const size_t numPrims = numPrimitives+4;
        const size_t numNodes = numPrimitives * BVH_NODE_PREALLOC_FACTOR + 4;
        bvh->init(numNodes,numPrims);
        node  = (BVHNode*) bvh->alloc_nodes->base();
        accel = (Triangle1*) bvh->alloc_tris->base();
        
        // === free previously allocated memory ===
        const size_t old_size_prims  = numPrimitivesOld * sizeof(PrimRef) + additional_size;
        if (prims) os_free(prims,old_size_prims);
        
        // === allocated memory for primrefs,nodes, and accel ===
        const size_t size_prims  = numPrimitives * sizeof(PrimRef) + additional_size;
        const size_t size_node  = numPrimitives * BVH_NODE_PREALLOC_FACTOR * sizeof(BVHNode) + additional_size;
        const size_t size_accel = numPrimitives * sizeof(Triangle1) + additional_size;
        numAllocatedNodes = size_node / sizeof(BVHNode);
        
        prims = (PrimRef      *) os_malloc(size_prims); memset(prims,0,size_prims);
        memset(node,0,size_node);
        memset(accel,0,size_accel);	
      }
    }
    
    void BVH4iBuilderFast::freeData()
    {
    }
    
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    
#if 1
    void BVH4iBuilderFast::computePrimRefs(const size_t threadID, const size_t numThreads)
    {
      const size_t numGroups = source->groups();
      const size_t startID = (threadID+0)*numPrimitives/numThreads;
      const size_t endID   = (threadID+1)*numPrimitives/numThreads;
      
      Scene* __restrict__ scene = (Scene*)geometry;
      PrimRef *__restrict__ const                 prims  = this->prims;
      
      // === find first group containing startID ===
      size_t g=0, numSkipped = 0;
      for (; g<numGroups; g++) {       
        Geometry* geom = scene->get(g);
        if (geom == NULL || geom->type != TRIANGLE_MESH) continue;
        TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
	if (unlikely(!mesh->isEnabled())) continue;
        const size_t numTriangles = mesh->numTriangles;
        if (numSkipped + numTriangles > startID) break;
        numSkipped += numTriangles;
      }
      
      // === start with first group containing startID ===
      Centroid_Scene_AABB bounds;
      bounds.reset();
      
      size_t num = 0;
      size_t currentID = startID;
      size_t offset = startID - numSkipped;
      for (; g<numGroups; g++) 
      {
        Geometry* geom = scene->get(g);
        if (geom == NULL || geom->type != TRIANGLE_MESH) continue;
        TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
	if (unlikely(!mesh->isEnabled())) continue;

        for (size_t i=offset; i<mesh->numTriangles && currentID < endID; i++, currentID++)	 
        { 			    
          const BBox3f bounds1 = mesh->bounds(i);
          bounds.extend(bounds1);
          prims[currentID] = PrimRef(bounds1,g,i);
        }
        if (currentID == endID) break;
        offset = 0;
      }
      global_bounds.extend_atomic(bounds);    
    }
    
#else
    
    void BVH4iBuilderFast::computePrimRefs(const size_t threadID, const size_t numThreads)
    {
      const size_t numGroups = source->groups();
      const size_t startID = (threadID+0)*numPrimitives/numThreads;
      const size_t endID   = (threadID+1)*numPrimitives/numThreads;
      
      PrimRef *__restrict__ const                 prims  = this->prims;
      
      // === find first group containing startID ===
      size_t g=0, numSkipped = 0;
      for (; g<numGroups; g++) {
        const size_t numTriangles = source->prims(g);
        if (numSkipped + numTriangles > startID) break;
        numSkipped += numTriangles;
      }
      
      // === start with first group containing startID ===
      Centroid_Scene_AABB bounds;
      bounds.reset();
      
      size_t num = 0;
      size_t currentID = startID;
      size_t offset = startID - numSkipped;
      for (; g<numGroups; g++) 
      {
        const size_t numTriangles = min(source->prims(g)-offset,endID-currentID);
        
        for (size_t i=0; i<numTriangles; i++, currentID++)	 
        { 			    
          const BBox3f bounds1 = source->bounds(g,offset+i);
          bounds.extend(bounds1);
          prims[currentID] = PrimRef(bounds1,g,offset+i);
        }
        if (currentID == endID) break;
        offset = 0;
      }
      global_bounds.extend_atomic(bounds);    
    }
#endif
    
    __forceinline void convertToSOALayoutBlock4(BVHNode  *bptr)
    {
      QBVHNode* qptr = (QBVHNode*)bptr;
      ssef min_x,min_y,min_z,min_a;
      ssef max_x,max_y,max_z,max_a;
      transpose((ssef)bptr[0].lower,(ssef)bptr[1].lower,(ssef)bptr[2].lower,(ssef)bptr[3].lower,min_x,min_y,min_z,min_a);
      
      const ssei min_d_org = cast(min_a);
      const sseb min_d_mask = bvhLeaf(min_d_org) != ssei(zero);
      const ssei min_d_node = qbvhCreateNode(bvhChildID(min_d_org)>>2,ssei(zero));
      const ssei min_d_leaf = (min_d_org ^ BVH_LEAF_MASK) | QBVH_LEAF_MASK;
      const ssei min_d = select(min_d_mask,min_d_leaf,min_d_node);
      
      transpose((ssef)bptr[0].upper,(ssef)bptr[1].upper,(ssef)bptr[2].upper,(ssef)bptr[3].upper,max_x,max_y,max_z,max_a);
      
      store4f_nt(qptr->min_x,min_x);
      store4f_nt(qptr->min_y,min_y);
      store4f_nt(qptr->min_z,min_z);
      store4i_nt(qptr->min_d,min_d);
      
      store4f_nt(qptr->max_x,max_x);
      store4f_nt(qptr->max_y,max_y);
      store4f_nt(qptr->max_z,max_z);
      store4i_nt(qptr->max_d,cast(max_a));       
    }
    
    void BVH4iBuilderFast::convertToSOALayout(const size_t threadID, const size_t numThreads)
    {
      const size_t startID = (threadID+0)*numNodes/numThreads;
      const size_t endID   = (threadID+1)*numNodes/numThreads;
      
      BVHNode  * __restrict__  bptr = ( BVHNode*)node + startID*4;
      
      for (size_t n=startID; n<endID; n++, bptr+=4)
        convertToSOALayoutBlock4(bptr);
      
    }
    
    __forceinline void computeAccelerationData(const size_t geomID,
                                               const size_t primID,     
                                               const Scene *__restrict__ const scene,
                                               Triangle1 * __restrict__ const acc)
    {
      const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
      const TriangleMeshScene::TriangleMesh::Triangle & tri = mesh->triangle(primID);
      
      const ssef v0 = select(0x7,(ssef) mesh->vertex(tri.v[0]),zero);
      const ssef v1 = select(0x7,(ssef) mesh->vertex(tri.v[1]),zero);
      const ssef v2 = select(0x7,(ssef) mesh->vertex(tri.v[2]),zero);
      
      const ssef e1 = v0 - v1;
      const ssef e2 = v2 - v0;	     
      const ssef normal = cross(e1,e2);
      
      store4f(&acc->v0,cast(insert<3>(cast(v0),primID)));
      store4f(&acc->v1,cast(insert<3>(cast(v1),geomID)));
      store4f(&acc->v2,cast(insert<3>(cast(v2),0)));
      store4f(&acc->Ng,cast(insert<3>(cast(normal),0)));
    }
    
    void BVH4iBuilderFast::createTriangle1(const size_t threadID, const size_t numThreads)
    {
      const size_t startID = (threadID+0)*numPrimitives/numThreads;
      const size_t endID   = (threadID+1)*numPrimitives/numThreads;
      
      Triangle1* __restrict__  acc = accel + startID;
      const PrimRef* __restrict__  bptr = prims + startID;
      
      for (size_t j=startID; j<endID; j++, bptr++, acc++)
        computeAccelerationData(bptr->geomID(),bptr->primID(),(Scene*)geometry,acc);
    }
    
    void BVH4iBuilderFast::buildSubTrees(const size_t threadID, const size_t numThreads)
    {
      while (true) 
      {
        BuildRecord br;
        if (!global_workStack.pop_largest(br)) // FIXME: might loose threads during build
        {
          /* global work queue empty => try to steal from neighboring queues */	  
          bool success = false;
          for (size_t i=0; i<numThreads; i++)
          {
            if (thread_workStack[(threadID+i)%numThreads].pop_smallest(br)) {
              success = true;
              break;
            }
          }
          /* found nothing to steal ? */
          if (!success) break; 
        }
        
        /* process local work queue */
	recurseSAH(br,RECURSE,threadID,numThreads);
        while (thread_workStack[threadID].pop_largest(br))
          recurseSAH(br,RECURSE,threadID,numThreads);
      }
    }
    
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    
    bool BVH4iBuilderFast::splitSequential(BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild)
    {
      /* mark as leaf if leaf threshold reached */
      if (current.items() <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD) {
        current.createLeaf();
        return false;
      }
      
      /* calculate binning function */
      Mapping<16> mapping(current.bounds);
      
      /* binning of centroids */
      Binner<16> binner;
      binner.bin(prims,current.begin,current.end,mapping);
      
      /* find best split */
      Split split; 
      binner.best(split,mapping);
      
      /* if we cannot find a valid split, enforce an arbitrary split */
      if (unlikely(split.pos == -1)) split_fallback(prims,current,leftChild,rightChild);
      
      /* partitioning of items */
      else binner.partition(prims, current.begin, current.end, split, mapping, leftChild, rightChild);
      
      if (leftChild.items()  <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD) leftChild.createLeaf();
      if (rightChild.items() <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD) rightChild.createLeaf();	
      return true;
    }
    
    bool BVH4iBuilderFast::splitParallel( BuildRecord &current,
                                          BuildRecord &leftChild,
                                          BuildRecord &rightChild,
                                          const size_t threadID,
                                          const size_t numThreads)
    {
      const unsigned int items = current.end - current.begin;
      assert(items >= BUILD_RECORD_SPLIT_THRESHOLD);
      
      /* mark as leaf if leaf threshold reached */
      if (items <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD) {
        current.createLeaf();
        return false;
      }
      
      /* parallel binning of centroids */
      parallelBinner.bin(current,prims,(PrimRef*)accel,threadID,numThreads);
      
      /* find best split */
      Split split; 
      parallelBinner.best(split);
      
      /* if we cannot find a valid split, enforce an arbitrary split */
      if (unlikely(split.pos == -1)) split_fallback(prims,current,leftChild,rightChild);
      
      /* parallel partitioning of items */
      else parallelBinner.partition((PrimRef*)accel,prims,split,leftChild,rightChild,threadID,numThreads);
      
      if (leftChild.items()  <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD) leftChild.createLeaf();
      if (rightChild.items() <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD) rightChild.createLeaf();
      return true;
    }
    
    __forceinline bool BVH4iBuilderFast::split(BuildRecord& current, BuildRecord& left, BuildRecord& right, const size_t mode, const size_t threadID, const size_t numThreads)
    {
      if (mode == BUILD_TOP_LEVEL && current.items() >= BUILD_RECORD_SPLIT_THRESHOLD)
        return splitParallel(current,left,right,threadID,numThreads);		  
      else
        return splitSequential(current,left,right);
    }
    
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    
    void BVH4iBuilderFast::createLeaf(BuildRecord& current, size_t threadIndex, size_t threadCount)
    {
#if defined(DEBUG)
      if (current.depth > BVH4i::maxBuildDepthLeaf) 
        throw std::runtime_error("ERROR: depth limit reached");
#endif
      
      /* create leaf for few primitives */
      if (current.items() <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD) {
        node[current.parentID].createLeaf(current.begin,current.items());
        return;
      }
      
      /* first split level */
      BuildRecord record0, record1;
      split_fallback(prims,current,record0,record1);
      
      /* second split level */
      BuildRecord children[4];
      split_fallback(prims,record0,children[0],children[1]);
      split_fallback(prims,record1,children[2],children[3]);
      
      /* allocate next four nodes */
      size_t numChildren = 4;
      const unsigned int currentIndex = allocNode(BVH4i::N);
      node[current.parentID].createNode(currentIndex,numChildren);
      
      /* recurse into each child */
      for (size_t i=0; i<numChildren; i++) 
      {
        node[currentIndex+i].lower = children[i].bounds.geometry.lower;
        node[currentIndex+i].upper = children[i].bounds.geometry.upper;
        children[i].parentID = currentIndex+i;
        children[i].depth = current.depth+1;
        createLeaf(children[i],threadIndex,threadCount);
      }
      
      //convertToSOALayoutBlock4((BVHNode*)node + currentIndex);
      
    }  
    
    __forceinline void BVH4iBuilderFast::recurse(BuildRecord& current, const size_t mode, const size_t threadID, const size_t numThreads)
    {
      if (mode == BUILD_TOP_LEVEL) {
        global_workStack.push_nolock(current);
      }
      else if (current.items() > THRESHOLD_FOR_SUBTREE_RECURSION) {
        if (!thread_workStack[threadID].push(current))
          recurseSAH(current,RECURSE,threadID,numThreads);
      }
      else
        recurseSAH(current,RECURSE,threadID,numThreads);
    }
    
    void BVH4iBuilderFast::recurseSAH(BuildRecord& current, const size_t mode, const size_t threadID, const size_t numThreads)
    {
      __aligned(64) BuildRecord children[BVH4i::N];
      
      /* create leaf node */
      if (current.depth >= BVH4i::maxBuildDepth || current.isLeaf()) {
        //node[current.parentID].createLeaf(current.begin,current.items());
        createLeaf(current,threadID,numThreads);
        return;
      }
      
      /* fill all 4 children by always splitting the one with the largest surface area */
      unsigned int numChildren = 1;
      children[0] = current;
      
      do {
        
        /* find best child with largest bounding box area */
        int bestChild = -1;
        float bestArea = neg_inf;
        for (unsigned int i=0; i<numChildren; i++)
        {
          /* ignore leaves as they cannot get split */
          if (children[i].isLeaf())
            continue;
          
          /* remember child with largest area */
          if (children[i].sceneArea() > bestArea) { 
            bestArea = children[i].sceneArea();
            bestChild = i;
          }
        }
        if (bestChild == -1) break;
        
        /*! split best child into left and right child */
        __aligned(64) BuildRecord left, right;
        if (!split(children[bestChild],left,right,mode,threadID,numThreads)) 
          continue;
        
        /* add new children left and right */
        left.depth = right.depth = current.depth+1;
        children[bestChild] = children[numChildren-1];
        children[numChildren-1] = left;
        children[numChildren+0] = right;
        numChildren++;
        
      } while (numChildren < BVH4i::N);
      
      /* create leaf node if no split is possible */
      if (numChildren == 1) {
        //node[current.parentID].createLeaf(current.begin,current.items());
        createLeaf(current,threadID,numThreads);
        return;
      }
      
      /* allocate next four nodes */
      const unsigned int currentIndex = allocNode(BVH4i::N);
      node[current.parentID].createNode(currentIndex,numChildren);
      
      /* recurse into each child */
      for (unsigned int i=0; i<numChildren; i++) 
      {
        node[currentIndex+i].lower = (Vec3fa) children[i].bounds.geometry.lower;
        node[currentIndex+i].upper = (Vec3fa) children[i].bounds.geometry.upper;
        children[i].parentID = currentIndex+i;
        recurse(children[i],mode,threadID,numThreads);
      }
      
      /* init unused nodes */
      const avxf init_node = load8f((float*)initQBVHNode);
      for (size_t i=numChildren; i<BVH4i::N; i++)
        store8f((float*)&node[currentIndex+i],init_node);
      
      //convertToSOALayoutBlock4((BVHNode*)node + currentIndex);
    }
    
    void BVH4iBuilderFast::build_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event) 
    {
      /* wait for all threads to enter */
      global_barrier.wait(threadIndex,threadCount);
      
      /* start measurement */
      double t0 = 0.0f;
      if (g_verbose >= 2) t0 = getSeconds();
      
      /* all worker threads enter tasking system */
      if (threadIndex != 0) {
        LockStepTaskScheduler::dispatchTaskMainLoop(threadIndex,threadCount); 
        return;
      }
      
      /* calculate list of primrefs */
      global_bounds.reset();
      LockStepTaskScheduler::dispatchTask( task_computePrimRefs, this, threadIndex, threadCount );
      
      /* allocate and initialize root node */
      atomicID.reset(BVH4i::N);
      node[0].lower = global_bounds.geometry.lower;
      node[0].upper = global_bounds.geometry.upper;
      
      /* create initial build record */
      BuildRecord br;
      br.init(global_bounds,0,numPrimitives);
      br.depth = 1;
      br.parentID = 0;
      
      /* initialize thread-local work stacks */
      for (size_t i=0;i<threadCount;i++)
        thread_workStack[i].reset();
      
      /* push initial build record to global work stack */
      global_workStack.reset();
      global_workStack.push_nolock(br);    
      
      /* work in multithreaded toplevel mode until sufficient subtasks got generated */
      while (global_workStack.size() < 4*threadCount && global_workStack.size()+BVH4i::N <= SIZE_WORK_STACK) 
      {
        BuildRecord br;
        if (!global_workStack.pop_nolock_largest(br)) break;
        recurseSAH(br,BUILD_TOP_LEVEL,threadIndex,threadCount);
      }
      
      /* now process all created subtasks on multiple threads */
      LockStepTaskScheduler::dispatchTask(task_buildSubTrees, this, threadIndex, threadCount );
      numNodes = atomicID >> 2;
      
      /* create triangle acceleration structure */
      LockStepTaskScheduler::dispatchTask( task_createTriangle1, this, threadIndex, threadCount );
      
      /* convert to SOA node layout */
      bvh->accel = accel;
      bvh->qbvh  = node;
      LockStepTaskScheduler::dispatchTask( task_convertToSOALayout, this, threadIndex, threadCount );
      
      //convertToSOALayoutBlock4((BVHNode*)node);
      
      /* update BVH4 */
      const QBVHNode* const __restrict__ qbvh  = (QBVHNode*)bvh->qbvh;
      bvh->root = qbvh[0].min_d[0]; 
      bvh->bounds = BBox3f(Vec3fa(qbvh->min_x[0],qbvh->min_y[0],qbvh->min_y[0]),
                           Vec3fa(qbvh->max_x[0],qbvh->max_y[0],qbvh->max_y[0]));
      
      /* release all threads again */
      LockStepTaskScheduler::releaseThreads(threadCount);
      
      /* stop measurement */
      if (g_verbose >= 2) dt = getSeconds()-t0;
    }

    Builder* BVH4iTriangle1BuilderObjectSplit4Fast (void* bvh, BuildSource* source, Scene* scene, const size_t minLeafSize, const size_t maxLeafSize) {
      return new BVH4iBuilderFast((BVH4i*)bvh,source,scene,minLeafSize,maxLeafSize);
    }
  }
}
