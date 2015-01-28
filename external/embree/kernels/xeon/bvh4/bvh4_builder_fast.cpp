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
#include "bvh4_builder_fast.h"
#include "bvh4_statistics.h"
#include "bvh4_builder_binner.h"

#include "geometry/triangle1.h"
#include "geometry/triangle4.h"
#include "geometry/triangle1v.h"
#include "geometry/triangle4v.h"

#define BVH_NODE_PREALLOC_FACTOR 1.0f
#define QBVH_BUILDER_LEAF_ITEM_THRESHOLD 4
#define THRESHOLD_FOR_SUBTREE_RECURSION 128
#define BUILD_RECORD_SPLIT_THRESHOLD 512

#define DBG(x) 

//#if defined(__USE_STAT_COUNTERS__)
//#define PROFILE
//#endif

namespace embree
{
  namespace isa
  {
    static double dt = 0.0f;

    std::auto_ptr<BVH4BuilderFast::GlobalState> BVH4BuilderFast::g_state(NULL);
    
    BVH4BuilderFast::BVH4BuilderFast (BVH4* bvh, BuildSource* source, Scene* scene, TriangleMeshScene::TriangleMesh* mesh, const size_t minLeafSize, const size_t maxLeafSize)
    : source(source), scene(scene), mesh(mesh), primTy(bvh->primTy), bvh(bvh), numGroups(0), numPrimitives(0), prims(NULL), bytesPrims(0), createSmallLeaf(NULL)
    {
      needAllThreads = true;
      if (mesh) needAllThreads = mesh->numTriangles > 50000;

      if      (&bvh->primTy == &SceneTriangle1::type ) createSmallLeaf = createTriangle1Leaf;
      else if (&bvh->primTy == &SceneTriangle4::type ) createSmallLeaf = createTriangle4Leaf;
      else if (&bvh->primTy == &SceneTriangle1v::type) createSmallLeaf = createTriangle1vLeaf;
      else if (&bvh->primTy == &SceneTriangle4v::type) createSmallLeaf = createTriangle4vLeaf;
      else if (&bvh->primTy == &TriangleMeshTriangle1::type ) createSmallLeaf = createTriangle1Leaf;
      else if (&bvh->primTy == &TriangleMeshTriangle4::type ) createSmallLeaf = createTriangle4Leaf;
      else if (&bvh->primTy == &TriangleMeshTriangle1v::type) createSmallLeaf = createTriangle1vLeaf;
      else if (&bvh->primTy == &TriangleMeshTriangle4v::type) createSmallLeaf = createTriangle4vLeaf;
      else throw std::runtime_error("BVH4BuilderFast: unknown primitive type");
    }
    
    BVH4BuilderFast::~BVH4BuilderFast () 
    {
      if (prims) os_free(prims,bytesPrims); prims = NULL;
      nodeAllocator.shrink(); 
      primAllocator.shrink();
      bvh->bytesNodes = nodeAllocator.bytesReserved;
      bvh->bytesPrimitives = primAllocator.bytesReserved;
    }
    
    void BVH4BuilderFast::build(size_t threadIndex, size_t threadCount) 
    {
      if (g_verbose >= 1)
        std::cout << "building BVH4 with " << TOSTRING(isa) "::BVH4BuilderFast ... " << std::flush;
      
      /* do some global inits first */
      init(threadIndex,threadCount);
      
#if defined(PROFILE)
      
      double dt_min = pos_inf;
      double dt_avg = 0.0f;
      double dt_max = neg_inf;
      for (size_t i=0; i<20; i++) 
      {
        if (!needAllThreads) {
          build_sequential(threadIndex,threadCount);
        } 
        else 
        {
          if (!g_state.get()) 
            g_state.reset(new GlobalState(threadCount));

          g_state->scheduler.init(threadCount);
          TaskScheduler::executeTask(threadIndex,threadCount,_build_parallel,this,threadCount,"build_parallel");
        }
        dt_min = min(dt_min,dt);
        dt_avg = dt_avg + dt;
        dt_max = max(dt_max,dt);
      }
      dt_avg /= double(20);
      
      std::cout << "[DONE]" << std::endl;
      std::cout << "  min = " << 1000.0f*dt_min << "ms (" << numPrimitives/dt_min*1E-6 << " Mtris/s)" << std::endl;
      std::cout << "  avg = " << 1000.0f*dt_avg << "ms (" << numPrimitives/dt_avg*1E-6 << " Mtris/s)" << std::endl;
      std::cout << "  max = " << 1000.0f*dt_max << "ms (" << numPrimitives/dt_max*1E-6 << " Mtris/s)" << std::endl;
      std::cout << "  node allocator = " 
                << 1E-6*nodeAllocator.next << " MB, " 
                << 1E-6*nodeAllocator.bytesAllocated << " MB, " 
                << 1E-6*nodeAllocator.bytesReserved << " MB" << std::endl;
      std::cout << "  primitive allocator = " 
                << 1E-6*primAllocator.next << " MB, " 
                << 1E-6*primAllocator.bytesAllocated << " MB, " 
                << 1E-6*primAllocator.bytesReserved << " MB" << std::endl;
      std::cout << BVH4Statistics(bvh).str();
      
#else
      
        if (!needAllThreads) {
          build_sequential(threadIndex,threadCount);
        } 
        else 
        {
          if (!g_state.get()) 
            g_state.reset(new GlobalState(threadCount));

          g_state->scheduler.init(threadCount);
          TaskScheduler::executeTask(threadIndex,threadCount,_build_parallel,this,threadCount,"build_parallel");
        }
      
      if (g_verbose >= 2) {
        double perf = numPrimitives/dt*1E-6;

        std::cout << "[DONE] " << 1000.0f*dt << "ms (" << perf << " Mtris/s)" << std::endl;
        std::cout << "  node allocator = " 
                  << 1E-6*nodeAllocator.next << " MB, " 
                  << 1E-6*nodeAllocator.bytesAllocated << " MB, " 
                  << 1E-6*nodeAllocator.bytesReserved << " MB" << std::endl;
        std::cout << "  primitive allocator = " 
                  << 1E-6*primAllocator.next << " MB, " 
                  << 1E-6*primAllocator.bytesAllocated << " MB, " 
                  << 1E-6*primAllocator.bytesReserved << " MB" << std::endl;
        std::cout << BVH4Statistics(bvh).str();
      }
#endif
    }
    
    void BVH4BuilderFast::init(size_t threadIndex, size_t threadCount)
    {
      //bvh->clear();
      bvh->init(0); // FIXME
      numGroups = scene->size();
      
      /* calculate size of scene */
      size_t numVertices = 0;
      size_t numPrimitivesOld = numPrimitives;
      if (mesh) {
        numPrimitives = mesh->numTriangles;
        numVertices = mesh->numVertices;
      }
      else {
        numPrimitives = 0;
        for (size_t i=0; i<scene->size(); i++) {
          Geometry* geom = scene->get(i);
          if (!geom || geom->type != TRIANGLE_MESH) continue;
          TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
          if (mesh->numTimeSteps != 1) continue;
          numPrimitives += mesh->numTriangles;
          numVertices += mesh->numVertices;
        }
      }
      bvh->numPrimitives = numPrimitives;
      if (primTy.needVertices) bvh->numVertices = numVertices;
      else                     bvh->numVertices = 0;
            
      if (numPrimitivesOld != numPrimitives)
      {
        /* free previously allocated memory */
        if (prims ) os_free(prims,bytesPrims);

        /* add one additional memory block for each thread in multi-threaded mode */
        size_t additionalBlocks = 1;
        if (needAllThreads) additionalBlocks = threadCount;
        
        /* allocate as much memory as likely needed and reserve conservative amounts of memory */
        size_t numPrimBlocks = bvh->primTy.blocks(numPrimitives);
        size_t numAllocatedNodes = min(size_t(0.6*numPrimBlocks),numPrimitives);
        size_t numAllocatedPrimitives = min(size_t(1.2*numPrimBlocks),numPrimitives);
#if defined(__X86_64__)
        size_t numReservedNodes = 2*numPrimitives;
        size_t numReservedPrimitives = 2*numPrimitives;
#else
        size_t numReservedNodes = 1.5*numAllocatedNodes;
        size_t numReservedPrimitives = 1.5*numAllocatedPrimitives;
#endif

        bytesPrims = numPrimitives * sizeof(PrimRef);
        size_t bytesAllocatedNodes      = numAllocatedNodes * sizeof(BVH4::Node);
        size_t bytesAllocatedPrimitives = numAllocatedPrimitives * bvh->primTy.bytes;
        bytesAllocatedPrimitives        = max(bytesAllocatedPrimitives,bytesPrims); // required as we store prims into primitive array for parallel splits
        size_t bytesReservedNodes       = numReservedNodes * sizeof(BVH4::Node);
        size_t bytesReservedPrimitives  = numReservedPrimitives * bvh->primTy.bytes;
        size_t blocksReservedNodes      = (bytesReservedNodes     +Allocator::blockSize-1)/Allocator::blockSize;
        size_t blocksReservedPrimitives = (bytesReservedPrimitives+Allocator::blockSize-1)/Allocator::blockSize;
        bytesReservedNodes      = Allocator::blockSize*(blocksReservedNodes      + additionalBlocks);
        bytesReservedPrimitives = Allocator::blockSize*(blocksReservedPrimitives + additionalBlocks);
        
        /* allocated memory for primrefs, nodes, and primitives */
        prims = (PrimRef* ) os_malloc(bytesPrims);  memset(prims,0,bytesPrims);
        nodeAllocator.init(bytesAllocatedNodes,bytesReservedNodes);
        primAllocator.init(bytesAllocatedPrimitives,bytesReservedPrimitives);
        
        bvh->nodes = nodeAllocator.data; 
        bvh->bytesNodes = nodeAllocator.bytesReserved;
        bvh->primitives = primAllocator.data; 
        bvh->bytesPrimitives = primAllocator.bytesReserved;
      }
    }
    
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    
#if 1
    void BVH4BuilderFast::computePrimRefs(const size_t threadID, const size_t numThreads)
    {
      const size_t startID = (threadID+0)*numPrimitives/numThreads;
      const size_t endID   = (threadID+1)*numPrimitives/numThreads;
      PrimRef *__restrict__ const                 prims  = this->prims;
      
      // === find first group containing startID ===
      size_t g=0, offset = 0;
      
      if (mesh) {
        g = mesh->id;
        offset = startID;
      }
      else
      {
        size_t numSkipped = 0;
        for (; g<numGroups; g++) {       
          Geometry* geom = scene->get(g);
          if (!geom || geom->type != TRIANGLE_MESH) continue;
          TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
          if (mesh->numTimeSteps != 1) continue;
          const size_t numTriangles = mesh->numTriangles;
          if (numSkipped + numTriangles > startID) break;
          numSkipped += numTriangles;
        }
        offset = startID - numSkipped;
      }
      
      // === start with first group containing startID ===
      Centroid_Scene_AABB bounds;
      bounds.reset();
      
      size_t num = 0;
      size_t currentID = startID;
      for (; g<numGroups; g++) 
      {
        Geometry* geom = scene->get(g);
        if (!geom || geom->type != TRIANGLE_MESH) continue;
        TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
        if (mesh->numTimeSteps != 1) continue;
        
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
    
    void BVH4BuilderFast::computePrimRefs(const size_t threadID, const size_t numThreads)
    {
      const size_t startID = (threadID+0)*numPrimitives/numThreads;
      const size_t endID   = (threadID+1)*numPrimitives/numThreads;
      
      PrimRef *__restrict__ const prims  = this->prims;
      
      // === find first group containing startID ===
      size_t g=0, numSkipped = 0;
      if (mesh) {
        g = mesh->id;
        offset = startID;
      }
      else 
      {
        for (; g<numGroups; g++) {
          const size_t numTriangles = source->prims(g);
          if (numSkipped + numTriangles > startID) break;
          numSkipped += numTriangles;
        }
        offset = startID - numSkipped;
      }
      
      // === start with first group containing startID ===
      Centroid_Scene_AABB bounds;
      bounds.reset();
      
      size_t num = 0;
      size_t currentID = startID;
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
    
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    
    void BVH4BuilderFast::createTriangle1Leaf(const BVH4BuilderFast* This, BuildRecord& current, Allocator& leafAlloc, size_t threadID)
    {
      size_t items = current.items();
      size_t start = current.begin;
      assert(items<=4);
      
      /* allocate leaf node */
      Triangle1* accel = (Triangle1*) leafAlloc.malloc(items*sizeof(Triangle1));
      *(NodeRef*)current.parentNode = This->bvh->encodeLeaf((char*)accel,items);
      
      for (size_t i=0; i<items; i++) 
      {	
        const size_t geomID = This->prims[start+i].geomID();
        const size_t primID = This->prims[start+i].primID();
        const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = This->scene->getTriangleMesh(geomID);
        const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
        
        const ssef v0 = select(0x7,(ssef)mesh->vertex(tri.v[0]),zero);
        const ssef v1 = select(0x7,(ssef)mesh->vertex(tri.v[1]),zero);
        const ssef v2 = select(0x7,(ssef)mesh->vertex(tri.v[2]),zero);
        
        const ssef e1 = v0 - v1;
        const ssef e2 = v2 - v0;	     
        const ssef normal = cross(e1,e2);
        
        store4f_nt(&accel[i].v0,cast(insert<3>(cast(v0),primID)));
        store4f_nt(&accel[i].v1,cast(insert<3>(cast(v1),geomID)));
        store4f_nt(&accel[i].v2,cast(insert<3>(cast(v2),mesh->mask)));
        store4f_nt(&accel[i].Ng,cast(insert<3>(cast(normal),0)));
      }
    }
    
    void BVH4BuilderFast::createTriangle4Leaf(const BVH4BuilderFast* This, BuildRecord& current, Allocator& leafAlloc, size_t threadID)
    {
      size_t items = current.items();
      size_t start = current.begin;
      assert(items<=4);
      
      /* allocate leaf node */
      Triangle4* accel = (Triangle4*) leafAlloc.malloc(sizeof(Triangle4));
      *(NodeRef*)current.parentNode = This->bvh->encodeLeaf((char*)accel,1);
      
      ssei vgeomID = -1, vprimID = -1, vmask = -1;
      sse3f v0 = zero, v1 = zero, v2 = zero;
      
      for (size_t i=0; i<items; i++)
      {
        const size_t geomID = This->prims[start+i].geomID();
        const size_t primID = This->prims[start+i].primID();
        const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = This->scene->getTriangleMesh(geomID);
        const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
        const Vec3fa& p0 = mesh->vertex(tri.v[0]);
        const Vec3fa& p1 = mesh->vertex(tri.v[1]);
        const Vec3fa& p2 = mesh->vertex(tri.v[2]);
        vgeomID [i] = geomID;
        vprimID [i] = primID;
        vmask   [i] = mesh->mask;
        v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
        v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
        v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
      }
      Triangle4::store_nt(accel,Triangle4(v0,v1,v2,vgeomID,vprimID,vmask));
    }
    
    void BVH4BuilderFast::createTriangle1vLeaf(const BVH4BuilderFast* This, BuildRecord& current, Allocator& leafAlloc, size_t threadID)
    {
      size_t items = current.items();
      size_t start = current.begin;
      assert(items<=4);
      
      /* allocate leaf node */
      Triangle1v* accel = (Triangle1v*) leafAlloc.malloc(items*sizeof(Triangle1v));
      *(NodeRef*)current.parentNode = This->bvh->encodeLeaf((char*)accel,items);
      
      for (size_t i=0; i<items; i++) 
      {	
        const size_t geomID = This->prims[start+i].geomID();
        const size_t primID = This->prims[start+i].primID();
        const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = This->scene->getTriangleMesh(geomID);
        const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
        
        const ssef v0 = select(0x7,(ssef)mesh->vertex(tri.v[0]),zero);
        const ssef v1 = select(0x7,(ssef)mesh->vertex(tri.v[1]),zero);
        const ssef v2 = select(0x7,(ssef)mesh->vertex(tri.v[2]),zero);
        
        const ssef e1 = v0 - v1;
        const ssef e2 = v2 - v0;	     
        const ssef normal = cross(e1,e2);
        
        store4f_nt(&accel[i].v0,cast(insert<3>(cast(v0),primID)));
        store4f_nt(&accel[i].v1,cast(insert<3>(cast(v1),geomID)));
        store4f_nt(&accel[i].v2,cast(insert<3>(cast(v2),mesh->mask)));
      }
    }
    
    void BVH4BuilderFast::createTriangle4vLeaf(const BVH4BuilderFast* This, BuildRecord& current, Allocator& leafAlloc, size_t threadID)
    {
      size_t items = current.items();
      size_t start = current.begin;
      assert(items<=4);
      
      /* allocate leaf node */
      Triangle4v* accel = (Triangle4v*) leafAlloc.malloc(sizeof(Triangle4v));
      *(NodeRef*)current.parentNode = This->bvh->encodeLeaf((char*)accel,1);
      
      ssei vgeomID = -1, vprimID = -1, vmask = -1;
      sse3f v0 = zero, v1 = zero, v2 = zero;
      
      for (size_t i=0; i<items; i++)
      {
        const size_t geomID = This->prims[start+i].geomID();
        const size_t primID = This->prims[start+i].primID();
        const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = This->scene->getTriangleMesh(geomID);
        const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
        const Vec3fa& p0 = mesh->vertex(tri.v[0]);
        const Vec3fa& p1 = mesh->vertex(tri.v[1]);
        const Vec3fa& p2 = mesh->vertex(tri.v[2]);
        vgeomID [i] = geomID;
        vprimID [i] = primID;
        vmask   [i] = mesh->mask;
        v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
        v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
        v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
      }
      Triangle4v::store_nt(accel,Triangle4v(v0,v1,v2,vgeomID,vprimID,vmask));
    }
    
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    
    bool BVH4BuilderFast::splitSequential(BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild)
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
    
    bool BVH4BuilderFast::splitParallel(BuildRecord &current, BuildRecord &leftChild, BuildRecord &rightChild, const size_t threadID, const size_t numThreads)
    {
      const unsigned int items = current.end - current.begin;
      assert(items >= BUILD_RECORD_SPLIT_THRESHOLD);
      
      /* mark as leaf if leaf threshold reached */
      if (items <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD) {
        current.createLeaf();
        return false;
      }
      
      /* use primitive array temporarily for parallel splits */
      PrimRef* tmp = (PrimRef*) primAllocator.base();

      /* parallel binning of centroids */
      g_state->parallelBinner.bin(current,prims,tmp,threadID,numThreads);
      
      /* find best split */
      Split split; 
      g_state->parallelBinner.best(split);
      
      /* if we cannot find a valid split, enforce an arbitrary split */
      if (unlikely(split.pos == -1)) split_fallback(prims,current,leftChild,rightChild);
      
      /* parallel partitioning of items */
      else g_state->parallelBinner.partition(tmp,prims,split,leftChild,rightChild,threadID,numThreads);
      
      if (leftChild.items()  <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD) leftChild.createLeaf();
      if (rightChild.items() <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD) rightChild.createLeaf();
      return true;
    }
    
    __forceinline bool BVH4BuilderFast::split(BuildRecord& current, BuildRecord& left, BuildRecord& right, const size_t mode, const size_t threadID, const size_t numThreads)
    {
      if (mode == BUILD_TOP_LEVEL && current.items() >= BUILD_RECORD_SPLIT_THRESHOLD)
        return splitParallel(current,left,right,threadID,numThreads);		  
      else
        return splitSequential(current,left,right);
    }
    
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    
    void BVH4BuilderFast::createLeaf(BuildRecord& current, Allocator& nodeAlloc, Allocator& leafAlloc, size_t threadIndex, size_t threadCount)
    {
#if defined(DEBUG)
      if (current.depth > BVH4::maxBuildDepthLeaf) 
        throw std::runtime_error("ERROR: depth limit reached");
#endif
      
      /* create leaf for few primitives */
      if (current.items() <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD) {
        createSmallLeaf(this,current,leafAlloc,threadIndex);
        return;
      }
      
      /* first split level */
      BuildRecord record0, record1;
      split_fallback(prims,current,record0,record1);
      
      /* second split level */
      BuildRecord children[4];
      split_fallback(prims,record0,children[0],children[1]);
      split_fallback(prims,record1,children[2],children[3]);

      /* allocate node */
      Node* node = (Node*) nodeAlloc.malloc(sizeof(Node)); node->clear();
      *(NodeRef*)current.parentNode = bvh->encodeNode(node);
      
      /* recurse into each child */
      for (size_t i=0; i<4; i++) 
      {
        node->set(i,children[i].bounds.geometry);
        children[i].parentNode = (size_t)&node->child(i);
        children[i].depth = current.depth+1;
        createLeaf(children[i],nodeAlloc,leafAlloc,threadIndex,threadCount);
      }
      BVH4::compact(node); // move empty nodes to the end
    }
    
    __forceinline void BVH4BuilderFast::recurse(BuildRecord& current, Allocator& nodeAlloc, Allocator& leafAlloc, const size_t mode, const size_t threadID, const size_t numThreads)
    {
      if (mode == BUILD_TOP_LEVEL) {
        g_state->workStack.push_nolock(current);
      }
      else if (mode == RECURSE_PARALLEL && current.items() > THRESHOLD_FOR_SUBTREE_RECURSION) {
        if (!g_state->threadStack[threadID].push(current))
          recurseSAH(current,nodeAlloc,leafAlloc,RECURSE_SEQUENTIAL,threadID,numThreads);
      }
      else
        recurseSAH(current,nodeAlloc,leafAlloc,mode,threadID,numThreads);
    }
    
    void BVH4BuilderFast::recurseSAH(BuildRecord& current, Allocator& nodeAlloc, Allocator& leafAlloc, const size_t mode, const size_t threadID, const size_t numThreads)
    {
      __aligned(64) BuildRecord children[BVH4::N];
      
      /* create leaf node */
      if (current.depth >= BVH4::maxBuildDepth || current.isLeaf()) {
        assert(mode != BUILD_TOP_LEVEL);
        createLeaf(current,nodeAlloc,leafAlloc,threadID,numThreads);
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
        
      } while (numChildren < BVH4::N);

      /* create leaf node if no split is possible */
      if (numChildren == 1) {
        assert(mode != BUILD_TOP_LEVEL);
        createLeaf(current,nodeAlloc,leafAlloc,threadID,numThreads);
        return;
      }
      
      /* allocate node */
      Node* node = (Node*) nodeAlloc.malloc(sizeof(Node)); node->clear();
      *(NodeRef*)current.parentNode = bvh->encodeNode(node);
      
      /* recurse into each child */
      for (unsigned int i=0; i<numChildren; i++) 
      {  
        node->set(i,children[i].bounds.geometry);
        children[i].parentNode = (size_t)&node->child(i);
        children[i].depth = current.depth+1;
        recurse(children[i],nodeAlloc,leafAlloc,mode,threadID,numThreads);
      }
    }
    
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    
    void BVH4BuilderFast::buildSubTrees(const size_t threadID, const size_t numThreads)
    {
      __aligned(64) Allocator nodeAlloc(nodeAllocator);
      __aligned(64) Allocator leafAlloc(primAllocator);
      
      while (true) 
      {
        BuildRecord br;
        if (!g_state->workStack.pop_largest(br)) // FIXME: might loose threads during build
        {
          /* global work queue empty => try to steal from neighboring queues */	  
          bool success = false;
          for (size_t i=0; i<numThreads; i++)
          {
            if (g_state->threadStack[(threadID+i)%numThreads].pop_smallest(br)) {
              success = true;
              break;
            }
          }
          /* found nothing to steal ? */
          if (!success) break; 
        }
        
        /* process local work queue */
	recurseSAH(br,nodeAlloc,leafAlloc,RECURSE_PARALLEL,threadID,numThreads);
        while (g_state->threadStack[threadID].pop_largest(br))
          recurseSAH(br,nodeAlloc,leafAlloc,RECURSE_PARALLEL,threadID,numThreads);
      }
    }

    void BVH4BuilderFast::build_sequential(size_t threadIndex, size_t threadCount) 
    {
      /* start measurement */
      double t0 = 0.0f;
      if (g_verbose >= 2) t0 = getSeconds();
      
      /* initialize node and leaf allocator */
      nodeAllocator.reset();
      primAllocator.reset();
      __aligned(64) Allocator nodeAlloc(nodeAllocator);
      __aligned(64) Allocator leafAlloc(primAllocator);
     
      /* create prim refs */
      global_bounds.reset();
      computePrimRefs(0,1);
      bvh->bounds = global_bounds.geometry;

      /* create initial build record */
      BuildRecord br;
      br.init(global_bounds,0,numPrimitives);
      br.depth = 1;
      br.parentNode = (size_t)&bvh->root;

      /* build BVH in single thread */
      recurseSAH(br,nodeAlloc,leafAlloc,RECURSE_SEQUENTIAL,threadIndex,threadCount);

      /* stop measurement */
      if (g_verbose >= 2) dt = getSeconds()-t0;
    }

    void BVH4BuilderFast::build_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event) 
    {
      /* wait for all threads to enter */
      g_state->barrier.wait(threadIndex,threadCount);
      
      /* start measurement */
      double t0 = 0.0f;
      if (g_verbose >= 2) t0 = getSeconds();
      
      /* all worker threads enter tasking system */
      if (threadIndex != 0) {
        g_state->scheduler.dispatchTaskMainLoop(threadIndex,threadCount); 
        return;
      }
      
      /* calculate list of primrefs */
      global_bounds.reset();
      g_state->scheduler.dispatchTask( task_computePrimRefs, this, threadIndex, threadCount );
      bvh->bounds = global_bounds.geometry;
      
      /* initialize node and leaf allocator */
      nodeAllocator.reset();
      primAllocator.reset();
      __aligned(64) Allocator nodeAlloc(nodeAllocator);
      __aligned(64) Allocator leafAlloc(primAllocator);

      /* create initial build record */
      BuildRecord br;
      br.init(global_bounds,0,numPrimitives);
      br.depth = 1;
      br.parentNode = (size_t)&bvh->root;
      
      /* initialize thread-local work stacks */
      for (size_t i=0; i<threadCount; i++)
        g_state->threadStack[i].reset();
      
      /* push initial build record to global work stack */
      g_state->workStack.reset();
      g_state->workStack.push_nolock(br);    
      
      /* work in multithreaded toplevel mode until sufficient subtasks got generated */
      while (g_state->workStack.size() < 4*threadCount && g_state->workStack.size()+BVH4::N <= SIZE_WORK_STACK) 
      {
        BuildRecord br;

        /* pop largest item for better load balancing */
        if (!g_state->workStack.pop_nolock_largest(br)) 
          break;
        
        /* guarantees to create no leaves in this stage */
        if (br.items() <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD)
          break;

        recurseSAH(br,nodeAlloc,leafAlloc,BUILD_TOP_LEVEL,threadIndex,threadCount);
      }
      
      /* now process all created subtasks on multiple threads */
      g_state->scheduler.dispatchTask(task_buildSubTrees, this, threadIndex, threadCount );
      
      /* release all threads again */
      g_state->scheduler.releaseThreads(threadCount);
      
      /* stop measurement */
      if (g_verbose >= 2) dt = getSeconds()-t0;
    }
    
    Builder* BVH4BuilderObjectSplit4Fast (void* bvh, BuildSource* source, Scene* scene, const size_t minLeafSize, const size_t maxLeafSize) {
      return new BVH4BuilderFast((BVH4*)bvh,source,scene,NULL,minLeafSize,maxLeafSize);
    }

    Builder* BVH4BuilderObjectSplit4TriangleMeshFast (void* bvh, TriangleMeshScene::TriangleMesh* mesh, const size_t minLeafSize, const size_t maxLeafSize) {
      return new BVH4BuilderFast((BVH4*)bvh,NULL,mesh->parent,mesh,minLeafSize,maxLeafSize);
    }
  }
}
