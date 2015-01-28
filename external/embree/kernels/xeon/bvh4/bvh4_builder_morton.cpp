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
#include "bvh4_builder_morton.h"
#include "bvh4_statistics.h"

#include "geometry/triangle1.h"
#include "geometry/triangle4.h"
#include "geometry/triangle1v.h"
#include "geometry/triangle4v.h"

#include "sys/tasklogger.h"

#define BVH_NODE_PREALLOC_FACTOR 1.1f

#define DBG(x) 

//#if defined(__USE_STAT_COUNTERS__)
//#define PROFILE
//#endif

namespace embree 
{
  namespace isa
  {
    static double dt = 0.0f;

    std::auto_ptr<BVH4BuilderMorton::MortonBuilderState> BVH4BuilderMorton::g_state(NULL);
    
    BVH4BuilderMorton::BVH4BuilderMorton (BVH4* bvh, BuildSource* source, Scene* scene, TriangleMeshScene::TriangleMesh* mesh, const size_t minLeafSize, const size_t maxLeafSize)
    : bvh(bvh), source(source), scene(scene), mesh(mesh), topLevelItemThreshold(0), encodeShift(0), encodeMask(0),
      morton(NULL), bytesMorton(0), numGroups(0), numPrimitives(0), numAllocatedPrimitives(0), numAllocatedNodes(0)
    {
      needAllThreads = true;
      if (mesh) needAllThreads = mesh->numTriangles > 50000;
      
      if (&bvh->primTy == &SceneTriangle1::type) {
        createSmallLeaf = createTriangle1Leaf;
        leafBounds = leafBoundsTriangle1;
      }
      else if (&bvh->primTy == &SceneTriangle4::type) {
        createSmallLeaf = createTriangle4Leaf;
        leafBounds = leafBoundsTriangle4;
      }
      else if (&bvh->primTy == &SceneTriangle1v::type) {
        createSmallLeaf = createTriangle1vLeaf;
        leafBounds = leafBoundsTriangle1v;
      }
      else if (&bvh->primTy == &SceneTriangle4v::type) {
        createSmallLeaf = createTriangle4vLeaf;
        leafBounds = leafBoundsTriangle4v;
      }
      else if (&bvh->primTy == &TriangleMeshTriangle1::type) {
        createSmallLeaf = createTriangle1Leaf;
        leafBounds = leafBoundsTriangle1;
      }
      else if (&bvh->primTy == &TriangleMeshTriangle4::type) {
        createSmallLeaf = createTriangle4Leaf;
        leafBounds = leafBoundsTriangle4;
      }
      else if (&bvh->primTy == &TriangleMeshTriangle1v::type) {
        createSmallLeaf = createTriangle1vLeaf;
        leafBounds = leafBoundsTriangle1v;
      }
      else if (&bvh->primTy == &TriangleMeshTriangle4v::type) {
        createSmallLeaf = createTriangle4vLeaf;
        leafBounds = leafBoundsTriangle4v;
      }
      else 
        throw std::runtime_error("BVH4BuilderMorton: unknown primitive type");
    }
    
    BVH4BuilderMorton::~BVH4BuilderMorton () 
    {
      if (morton) os_free(morton,bytesMorton);
      nodeAllocator.shrink(); 
      primAllocator.shrink();
      bvh->bytesNodes = nodeAllocator.bytesReserved;
      bvh->bytesPrimitives = primAllocator.bytesReserved;
    }
    
    void BVH4BuilderMorton::build(size_t threadIndex, size_t threadCount) 
    {
      if (g_verbose >= 2)
        std::cout << "building BVH4 with " << TOSTRING(isa) << "::BVH4BuilderMorton ... " << std::flush;
      
      /* do some global inits first */
      init(threadIndex,threadCount);
      
#if defined(PROFILE)
      
      double dt_min = pos_inf;
      double dt_avg = 0.0f;
      double dt_max = neg_inf;
      for (size_t i=0; i<200; i++) 
      {
        if (needAllThreads) 
        {
          if (!g_state.get()) g_state.reset(new MortonBuilderState);
          scheduler.init(threadCount);
          TaskScheduler::executeTask(threadIndex,threadCount,_build_parallel_morton,this,threadCount,"build_parallel_morton");
        } else {
          build_sequential_morton(threadIndex,threadCount);
        }
        dt_min = min(dt_min,dt);
        dt_avg = dt_avg + dt;
        dt_max = max(dt_max,dt);
      }
      dt_avg /= double(200);
      
      std::cout << "[DONE]" << std::endl;
      std::cout << "  min = " << 1000.0f*dt_min << "ms (" << numPrimitives/dt_min*1E-6 << " Mtris/s)" << std::endl;
      std::cout << "  avg = " << 1000.0f*dt_avg << "ms (" << numPrimitives/dt_avg*1E-6 << " Mtris/s)" << std::endl;
      std::cout << "  max = " << 1000.0f*dt_max << "ms (" << numPrimitives/dt_max*1E-6 << " Mtris/s)" << std::endl;
      std::cout << BVH4Statistics(bvh).str();
      
#else
      
      if (needAllThreads) 
      {
        if (!g_state.get()) g_state.reset(new MortonBuilderState);
        scheduler.init(threadCount);
        TaskScheduler::executeTask(threadIndex,threadCount,_build_parallel_morton,this,threadCount,"build_parallel_morton");
      } else {
        build_sequential_morton(threadIndex,threadCount);
      }

      if (g_verbose >= 2) {
        double perf = numPrimitives/dt*1E-6;
        std::cout << "[DONE] " << 1000.0f*dt << "ms (" << perf << " Mtris/s)" << std::endl;
        std::cout << BVH4Statistics(bvh).str();
      }
#endif
    }
    
    void BVH4BuilderMorton::init(size_t threadIndex, size_t threadCount)
    {
      //bvh->clear();
      bvh->init(numPrimitives);
      
      /* calculate size of scene */
      size_t numVertices = 0;
      size_t numPrimitivesOld = numPrimitives;
      numGroups = scene->size();
      if (mesh) {
        numPrimitives = mesh->numTriangles;
        numVertices = mesh->numVertices;
      }
      else {
        numPrimitives = 0;
        for (size_t i=0; i<numGroups; i++) 
        {
          Geometry* geom = scene->get(i);
          if (!geom || geom->type != TRIANGLE_MESH) continue;
          TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
          if (mesh->numTimeSteps != 1) continue;
          numPrimitives += mesh->numTriangles;
          numVertices += mesh->numVertices;
        }
      }
      bvh->numPrimitives = numPrimitives;
      if (bvh->primTy.needVertices) bvh->numVertices = numVertices;
      else                          bvh->numVertices = 0;
      
      size_t maxPrimsPerGroup = 0;
      if (mesh) 
        maxPrimsPerGroup = numPrimitives;
      else {
        for (size_t g=0; g<numGroups; g++) 
          maxPrimsPerGroup = max(maxPrimsPerGroup,source->prims(g));
      }

      /* calculate groupID, primID encoding */
      encodeShift = __bsr(maxPrimsPerGroup) + 1;
      encodeMask = ((size_t)1 << encodeShift)-1;
      size_t maxGroups = ((size_t)1 << (31-encodeShift))-1;
      
      if (maxPrimsPerGroup > encodeMask || numGroups > maxGroups)
      {
        std::cout << "ENCODING ERROR" << std::endl;
        DBG_PRINT(numGroups);
        DBG_PRINT(numPrimitives);
        DBG_PRINT(maxPrimsPerGroup);
        DBG_PRINT(encodeMask);
        DBG_PRINT(maxGroups);
        exit(1);
      }
      
      if (mesh) {
        numPrimitives = mesh->numTriangles;
      }
      
      /* preallocate arrays */
      if (numPrimitivesOld != numPrimitives)
      {
        /* free previously allocated memory */
        if (morton) os_free(morton,bytesMorton);

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
        bytesMorton = ((numPrimitives+7)&(-8)) * sizeof(MortonID32Bit);
        size_t bytesAllocatedNodes      = numAllocatedNodes * sizeof(BVH4::Node);
        size_t bytesAllocatedPrimitives = numAllocatedPrimitives * bvh->primTy.bytes;
        size_t bytesReservedNodes       = numReservedNodes * sizeof(BVH4::Node);
        size_t bytesReservedPrimitives  = numReservedPrimitives * bvh->primTy.bytes;
        size_t blocksReservedNodes      = (bytesReservedNodes     +Allocator::blockSize-1)/Allocator::blockSize;
        size_t blocksReservedPrimitives = (bytesReservedPrimitives+Allocator::blockSize-1)/Allocator::blockSize;
        bytesReservedNodes      = Allocator::blockSize*(blocksReservedNodes      + additionalBlocks);
        bytesReservedPrimitives = Allocator::blockSize*(blocksReservedPrimitives + additionalBlocks);
        bytesAllocatedNodes = max(bytesAllocatedNodes,bytesMorton); 
        bytesReservedNodes  = max(bytesReservedNodes,bytesMorton); 

        /* allocated memory for primrefs, nodes, and primitives */
        morton = (MortonID32Bit* ) os_malloc(bytesMorton); memset(morton,0,bytesMorton);
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
    
    void BVH4BuilderMorton::initThreadState(const size_t threadID, const size_t numThreads)
    {
      const size_t startID = (threadID+0)*numPrimitives/numThreads;
      const size_t endID   = (threadID+1)*numPrimitives/numThreads;
      
      if (mesh) 
      {
        /* store start group and offset */
        g_state->startGroup[threadID] = mesh->id;
        g_state->startGroupOffset[threadID] = startID;
      }
      else
      {
        /* find first group containing startID */
        size_t group = 0, skipped = 0;
        for (; group<numGroups; group++) 
        {       
          Geometry* geom = scene->get(group);
          if (!geom || geom->type != TRIANGLE_MESH) continue;
          TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
          if (mesh->numTimeSteps != 1) continue;
          const size_t numTriangles = mesh->numTriangles;
          if (skipped + numTriangles > startID) break; 
         skipped += numTriangles;
        }
        
        /* store start group and offset */
        g_state->startGroup[threadID] = group;
        g_state->startGroupOffset[threadID] = startID - skipped;
      }
    }

    Centroid_Scene_AABB BVH4BuilderMorton::computeBounds()
    {
      Centroid_Scene_AABB bounds;
      bounds.reset();

      if (mesh) 
      {
        for (size_t i=0; i<mesh->numTriangles; i++)	 
          bounds.extend(mesh->bounds(i));
      }
      else
      {
        for (size_t group=0; group<numGroups; group++) 
        {       
          Geometry* geom = scene->get(group);
          if (!geom || geom->type != TRIANGLE_MESH) continue;
          TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
          if (mesh->numTimeSteps != 1) continue;
          for (size_t i=0; i<mesh->numTriangles; i++)	 
            bounds.extend(mesh->bounds(i));
        }
      }
      return bounds;
    }
    
    void BVH4BuilderMorton::computeBounds(const size_t threadID, const size_t numThreads)
    {
      const size_t startID = (threadID+0)*numPrimitives/numThreads;
      const size_t endID   = (threadID+1)*numPrimitives/numThreads;
      
      Centroid_Scene_AABB bounds;
      bounds.reset();
      
      size_t currentID = startID;
      size_t offset = g_state->startGroupOffset[threadID];
      
      for (size_t group = g_state->startGroup[threadID]; group<numGroups; group++) 
      {       
        Geometry* geom = scene->get(group);
        if (!geom || geom->type != TRIANGLE_MESH) continue;
        TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
        if (mesh->numTimeSteps != 1) continue;
        
        for (size_t i=offset; i<mesh->numTriangles && currentID < endID; i++, currentID++)	 
          bounds.extend(mesh->bounds(i));
        
        offset = 0;
        if (currentID == endID) break;
      }
      
      global_bounds.extend_atomic(bounds);    
    }

    void BVH4BuilderMorton::computeMortonCodes(const size_t startID, const size_t endID, 
                                               const size_t startGroup, const size_t startOffset, 
                                               MortonID32Bit* __restrict__ const dest)
    {
      /* compute mapping from world space into 3D grid */
      const ssef base     = (ssef)global_bounds.centroid2.lower;
      const ssef diagonal = (ssef)global_bounds.centroid2.upper - (ssef)global_bounds.centroid2.lower;
      const ssef scale    = select(diagonal != 0, rcp(diagonal) * ssef(LATTICE_SIZE_PER_DIM * 0.99f),ssef(0.0f));
      
      size_t currentID = startID;
      size_t offset = startOffset;
      
      /* use SSE to calculate morton codes */
      size_t slots = 0;
      ssei ax = 0, ay = 0, az = 0, ai = 0;
      
      for (size_t group = startGroup; group<numGroups; group++) 
      {       
        Geometry* geom = scene->get(group);
        if (!geom || geom->type != TRIANGLE_MESH) continue;
        TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
        if (mesh->numTimeSteps != 1) continue;
        const size_t numTriangles = min(mesh->numTriangles-offset,endID-currentID);
        
        for (size_t i=0; i<numTriangles; i++)	  
        {
          const BBox3f b = mesh->bounds(offset+i);
          const ssef lower = (ssef)b.lower;
          const ssef upper = (ssef)b.upper;
          const ssef centroid = lower+upper;
          const ssei binID = ssei((centroid-base)*scale);
          const unsigned int index = (group << encodeShift) | (offset+i);
          ax[slots] = extract<0>(binID);
          ay[slots] = extract<1>(binID);
          az[slots] = extract<2>(binID);
          ai[slots] = index;
          slots++;
          currentID++;
          
          if (slots == 4)
          {
            const ssei code = bitInterleave(ax,ay,az);
            storeu4i(&dest[currentID-4],unpacklo(code,ai));
            storeu4i(&dest[currentID-2],unpackhi(code,ai));
            slots = 0;
          }
        }
        offset = 0;
        if (currentID == endID) break;
      }
      
      if (slots != 0)
      {
        const ssei code = bitInterleave(ax,ay,az);
        for (size_t i=0; i<slots; i++) {
          dest[currentID-slots+i].index = ai[i];
          dest[currentID-slots+i].code = code[i];
        }
      }
    }
    
    void BVH4BuilderMorton::computeMortonCodes(const size_t threadID, const size_t numThreads)
    {      
      const size_t startID = (threadID+0)*numPrimitives/numThreads;
      const size_t endID   = (threadID+1)*numPrimitives/numThreads;
      
      /* store the morton codes temporarily in 'node' memory */
      MortonID32Bit* __restrict__ const dest = (MortonID32Bit*)nodeAllocator.data; 
      computeMortonCodes(startID,endID,g_state->startGroup[threadID],g_state->startGroupOffset[threadID],dest);
    }
    
    void BVH4BuilderMorton::recreateMortonCodes(SmallBuildRecord& current) const
    {
      assert(current.size() > 4);
      Centroid_Scene_AABB global_bounds;
      global_bounds.reset();
      
      for (size_t i=current.begin; i<current.end; i++)
      {
        const size_t index  = morton[i].index;
        const size_t primID = index & encodeMask; 
        const size_t geomID = index >> encodeShift; 
        global_bounds.extend(scene->getTriangleMesh(geomID)->bounds(primID));
      }
      
      /* compute mapping from world space into 3D grid */
      const ssef base     = (ssef)global_bounds.centroid2.lower;
      const ssef diagonal = (ssef)global_bounds.centroid2.upper - (ssef)global_bounds.centroid2.lower;
      const ssef scale    = select(diagonal != 0,rcp(diagonal) * ssef(LATTICE_SIZE_PER_DIM * 0.99f),ssef(0.0f));
      
      for (size_t i=current.begin; i<current.end; i++)
      {
        const size_t index  = morton[i].index;
        const size_t primID = index & encodeMask; 
        const size_t geomID = index >> encodeShift; 
        const BBox3f b = scene->getTriangleMesh(geomID)->bounds(primID);
        const ssef lower = (ssef)b.lower;
        const ssef upper = (ssef)b.upper;
        const ssef centroid = lower+upper;
        const ssei binID = ssei((centroid-base)*scale);
        const unsigned int bx = extract<0>(binID);
        const unsigned int by = extract<1>(binID);
        const unsigned int bz = extract<2>(binID);
        const unsigned int code = bitInterleave(bx,by,bz);
        morton[i].code  = code;
      }
      quicksort_insertionsort_ascending<MortonID32Bit,512>(morton,current.begin,current.end-1); 
      
#if defined(DEBUG)
      for (size_t i=current.begin; i<current.end-1; i++)
        assert(morton[i].code <= morton[i+1].code);
#endif	    
    }
    
    void BVH4BuilderMorton::radixsort(const size_t threadID, const size_t numThreads)
    {
      //size_t taskID = TaskLogger::beginTask(threadID,"BVH4BuilderMorton::radixsort",0);

      const size_t startID = (threadID+0)*numPrimitives/numThreads;
      const size_t endID   = (threadID+1)*numPrimitives/numThreads;
      
      MortonID32Bit* __restrict__ mortonID[2];
      mortonID[0] = (MortonID32Bit*) morton; 
      //mortonID[1] = (MortonID32Bit*) node;
      mortonID[1] = (MortonID32Bit*) nodeAllocator.data;
      MortonBuilderState::ThreadRadixCountTy* radixCount = g_state->radixCount;
      
      /* we need 3 iterations to process all 32 bits */
      for (size_t b=0; b<3; b++)
      {
        const MortonID32Bit* __restrict src = (MortonID32Bit*) &mortonID[((b+1)%2)][0];
        MortonID32Bit*       __restrict dst = (MortonID32Bit*) &mortonID[((b+0)%2)][0];
        
        /* shift and mask to extract some number of bits */
        const unsigned int mask = RADIX_BUCKETS_MASK;
        const unsigned int shift = b * RADIX_BITS;
        
        /* count how many items go into the buckets */
        for (size_t i=0; i<RADIX_BUCKETS; i++)
          radixCount[threadID][i] = 0;
        
        for (size_t i=startID; i<endID; i++) {
          const size_t index = src[i].get(shift, mask);
          radixCount[threadID][index]++;
        }
        scheduler.syncThreads(threadID,numThreads);
        
        /* calculate total number of items for each bucket */
        __aligned(64) size_t total[RADIX_BUCKETS];
        for (size_t i=0; i<RADIX_BUCKETS; i++)
          total[i] = 0;
        
        for (size_t i=0; i<numThreads; i++)
          for (size_t j=0; j<RADIX_BUCKETS; j++)
            total[j] += radixCount[i][j];
        
        /* calculate start offset of each bucket */
        __aligned(64) size_t offset[RADIX_BUCKETS];
        offset[0] = 0;
        for (size_t i=1; i<RADIX_BUCKETS; i++)    
          offset[i] = offset[i-1] + total[i-1];
        
        /* calculate start offset of each bucket for this thread */
        for (size_t j=0; j<RADIX_BUCKETS; j++)
          for (size_t i=0; i<threadID; i++)
            offset[j] += radixCount[i][j];
        
        /* copy items into their buckets */
        for (size_t i=startID; i<endID; i++) {
          const size_t index = src[i].get(shift, mask);
          dst[offset[index]++] = src[i];
        }
        if (b < 2) scheduler.syncThreads(threadID,numThreads);
      }

      //TaskLogger::endTask(threadID,taskID);
    }
    
    void BVH4BuilderMorton::recurseSubMortonTrees(const size_t threadID, const size_t numThreads)
    {
      __aligned(64) Allocator nodeAlloc(nodeAllocator);
      __aligned(64) Allocator leafAlloc(primAllocator);
      while (true)
      {
        const unsigned int taskID = scheduler.taskCounter.inc();
        if (taskID >= g_state->numBuildRecords) break;
        
        //size_t id = TaskLogger::beginTask(threadID,"BVH4BuilderMorton::subtree",0);
        recurse(g_state->buildRecords[taskID],nodeAlloc,leafAlloc,RECURSE,threadID);
        g_state->buildRecords[taskID].parent->setBarrier();
        g_state->workStack.push(g_state->buildRecords[taskID]);
        //TaskLogger::endTask(threadID,id);
      }
    }
    
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    
    void BVH4BuilderMorton::createTriangle1Leaf(const BVH4BuilderMorton* This, SmallBuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3f& box_o)
    {
      ssef lower(pos_inf);
      ssef upper(neg_inf);
      size_t items = current.size();
      size_t start = current.begin;
      assert(items<=4);
      
      /* allocate leaf node */
      Triangle1* accel = (Triangle1*) leafAlloc.malloc(items*sizeof(Triangle1));
      *current.parent = This->bvh->encodeLeaf((char*)accel,items);
      
      for (size_t i=0; i<items; i++) 
      {	
        const size_t index = This->morton[start+i].index;
        const size_t primID = index & This->encodeMask; 
        const size_t geomID = index >> This->encodeShift; 
        const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = This->scene->getTriangleMesh(geomID);
        const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
        
        const ssef v0 = select(0x7,(ssef)mesh->vertex(tri.v[0]),zero);
        const ssef v1 = select(0x7,(ssef)mesh->vertex(tri.v[1]),zero);
        const ssef v2 = select(0x7,(ssef)mesh->vertex(tri.v[2]),zero);
        
        lower = min(lower,v0,v1,v2);
        upper = max(upper,v0,v1,v2);
        
        const ssef e1 = v0 - v1;
        const ssef e2 = v2 - v0;	     
        const ssef normal = cross(e1,e2);
        
        store4f_nt(&accel[i].v0,cast(insert<3>(cast(v0),primID)));
        store4f_nt(&accel[i].v1,cast(insert<3>(cast(v1),geomID)));
        store4f_nt(&accel[i].v2,cast(insert<3>(cast(v2),mesh->mask)));
        store4f_nt(&accel[i].Ng,cast(insert<3>(cast(normal),0)));
      }
      box_o = BBox3f((Vec3fa)lower,(Vec3fa)upper);
    }
    
    void BVH4BuilderMorton::createTriangle4Leaf(const BVH4BuilderMorton* This, SmallBuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3f& box_o)
    {
      ssef lower(pos_inf);
      ssef upper(neg_inf);
      size_t items = current.size();
      size_t start = current.begin;
      assert(items<=4);
      
      /* allocate leaf node */
      Triangle4* accel = (Triangle4*) leafAlloc.malloc(sizeof(Triangle4));
      *current.parent = This->bvh->encodeLeaf((char*)accel,1);
      
      ssei vgeomID = -1, vprimID = -1, vmask = -1;
      sse3f v0 = zero, v1 = zero, v2 = zero;
      
      for (size_t i=0; i<items; i++)
      {
        const size_t index = This->morton[start+i].index;
        const size_t primID = index & This->encodeMask; 
        const size_t geomID = index >> This->encodeShift; 
        const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = This->scene->getTriangleMesh(geomID);
        const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
        const Vec3fa& p0 = mesh->vertex(tri.v[0]);
        const Vec3fa& p1 = mesh->vertex(tri.v[1]);
        const Vec3fa& p2 = mesh->vertex(tri.v[2]);
        lower = min(lower,(ssef)p0,(ssef)p1,(ssef)p2);
        upper = max(upper,(ssef)p0,(ssef)p1,(ssef)p2);
        vgeomID [i] = geomID;
        vprimID [i] = primID;
        vmask   [i] = mesh->mask;
        v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
        v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
        v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
      }
      Triangle4::store_nt(accel,Triangle4(v0,v1,v2,vgeomID,vprimID,vmask));
      //new (accel) Triangle4(v0,v1,v2,vgeomID,vprimID,vmask);
      box_o = BBox3f((Vec3fa)lower,(Vec3fa)upper);
    }
    
    void BVH4BuilderMorton::createTriangle1vLeaf(const BVH4BuilderMorton* This, SmallBuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3f& box_o)
    {
      ssef lower(pos_inf);
      ssef upper(neg_inf);
      size_t items = current.size();
      size_t start = current.begin;
      assert(items<=4);
      
      /* allocate leaf node */
      Triangle1v* accel = (Triangle1v*) leafAlloc.malloc(items*sizeof(Triangle1v));
      *current.parent = This->bvh->encodeLeaf((char*)accel,items);
      
      for (size_t i=0; i<items; i++) 
      {	
        const size_t index = This->morton[start+i].index;
        const size_t primID = index & This->encodeMask; 
        const size_t geomID = index >> This->encodeShift; 
        const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = This->scene->getTriangleMesh(geomID);
        const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
        
        const ssef v0 = select(0x7,(ssef)mesh->vertex(tri.v[0]),zero);
        const ssef v1 = select(0x7,(ssef)mesh->vertex(tri.v[1]),zero);
        const ssef v2 = select(0x7,(ssef)mesh->vertex(tri.v[2]),zero);
        
        lower = min(lower,v0,v1,v2);
        upper = max(upper,v0,v1,v2);
        
        const ssef e1 = v0 - v1;
        const ssef e2 = v2 - v0;	     
        const ssef normal = cross(e1,e2);
        
        store4f_nt(&accel[i].v0,cast(insert<3>(cast(v0),primID)));
        store4f_nt(&accel[i].v1,cast(insert<3>(cast(v1),geomID)));
        store4f_nt(&accel[i].v2,cast(insert<3>(cast(v2),mesh->mask)));
      }
      box_o = BBox3f((Vec3fa)lower,(Vec3fa)upper);
    }
    
    void BVH4BuilderMorton::createTriangle4vLeaf(const BVH4BuilderMorton* This, SmallBuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3f& box_o)
    {
      ssef lower(pos_inf);
      ssef upper(neg_inf);
      size_t items = current.size();
      size_t start = current.begin;
      assert(items<=4);
      
      /* allocate leaf node */
      Triangle4v* accel = (Triangle4v*) leafAlloc.malloc(sizeof(Triangle4v));
      *current.parent = This->bvh->encodeLeaf((char*)accel,1);
      
      ssei vgeomID = -1, vprimID = -1, vmask = -1;
      sse3f v0 = zero, v1 = zero, v2 = zero;
      
      for (size_t i=0; i<items; i++)
      {
        const size_t index = This->morton[start+i].index;
        const size_t primID = index & This->encodeMask; 
        const size_t geomID = index >> This->encodeShift; 
        const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = This->scene->getTriangleMesh(geomID);
        const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
        const Vec3fa& p0 = mesh->vertex(tri.v[0]);
        const Vec3fa& p1 = mesh->vertex(tri.v[1]);
        const Vec3fa& p2 = mesh->vertex(tri.v[2]);
        lower = min(lower,(ssef)p0,(ssef)p1,(ssef)p2);
        upper = max(upper,(ssef)p0,(ssef)p1,(ssef)p2);
        vgeomID [i] = geomID;
        vprimID [i] = primID;
        vmask   [i] = mesh->mask;
        v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
        v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
        v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
      }
      Triangle4v::store_nt(accel,Triangle4v(v0,v1,v2,vgeomID,vprimID,vmask));
      //new (accel) Triangle4v(v0,v1,v2,vgeomID,vprimID,vmask);
      box_o = BBox3f((Vec3fa)lower,(Vec3fa)upper);
    }
    
    void BVH4BuilderMorton::split_fallback(SmallBuildRecord& current, SmallBuildRecord& leftChild, SmallBuildRecord& rightChild) const
    {
      const unsigned int center = (current.begin + current.end)/2;
      leftChild.init(current.begin,center);
      rightChild.init(center,current.end);
    }
    
    BBox3f BVH4BuilderMorton::createLeaf(SmallBuildRecord& current, Allocator& nodeAlloc, Allocator& leafAlloc, size_t threadID)
    {
#if defined(DEBUG)
      if (current.depth > BVH4::maxBuildDepthLeaf) 
        throw std::runtime_error("ERROR: depth limit reached");
#endif
      
      /* create leaf for few primitives */
      if (current.size() <= MORTON_LEAF_THRESHOLD) {
        BBox3f bounds;
        createSmallLeaf(this,current,leafAlloc,threadID,bounds);
        return bounds;
      }
      
      /* first split level */
      SmallBuildRecord record0, record1;
      split_fallback(current,record0,record1);
      
      /* second split level */
      SmallBuildRecord children[4];
      split_fallback(record0,children[0],children[1]);
      split_fallback(record1,children[2],children[3]);
      
      /* allocate node */
      Node* node = (Node*) nodeAlloc.malloc(sizeof(Node)); node->clear();
      *current.parent = bvh->encodeNode(node);
      
      /* recurse into each child */
      BBox3f bounds0 = empty;
      for (size_t i=0; i<4; i++) {
        children[i].parent = &node->child(i);
        children[i].depth = current.depth+1;
        BBox3f bounds = createLeaf(children[i],nodeAlloc,leafAlloc,threadID);
        bounds0.extend(bounds);
        node->set(i,bounds);
      }
      BVH4::compact(node); // move empty nodes to the end
      return bounds0;
    }  
    
    __forceinline void BVH4BuilderMorton::split(SmallBuildRecord& current,
                                                SmallBuildRecord& left,
                                                SmallBuildRecord& right) const
    {
      const unsigned int code_start = morton[current.begin].code;
      const unsigned int code_end   = morton[current.end-1].code;
      unsigned int bitpos = clz(code_start^code_end);
      
      /* if all items mapped to same morton code, then create new morton codes for the items */
      if (unlikely(bitpos == 32)) 
      {
        recreateMortonCodes(current);
        const unsigned int code_start = morton[current.begin].code;
        const unsigned int code_end   = morton[current.end-1].code;
        bitpos = clz(code_start^code_end);
        
        /* if the morton code is still the same, goto fall back split */
        if (unlikely(bitpos == 32)) 
        {
          size_t center = (current.begin + current.end)/2; 
          left.init(current.begin,center);
          right.init(center,current.end);
          return;
        }
      }
      
      /* split the items at the topmost different morton code bit */
      const unsigned int bitpos_diff = 31-bitpos;
      const unsigned int bitmask = 1 << bitpos_diff;
      
      /* find location where bit differs using binary search */
      size_t begin = current.begin;
      size_t end   = current.end;
      while (begin + 1 != end) {
        const size_t mid = (begin+end)/2;
        const unsigned bit = morton[mid].code & bitmask;
        if (bit == 0) begin = mid; else end = mid;
      }
      size_t center = end;
#if defined(DEBUG)      
      for (unsigned int i=begin;  i<center; i++) assert((morton[i].code & bitmask) == 0);
      for (unsigned int i=center; i<end;    i++) assert((morton[i].code & bitmask) == bitmask);
#endif
      
      left.init(current.begin,center);
      right.init(center,current.end);
    }
    
    BBox3f BVH4BuilderMorton::recurse(SmallBuildRecord& current, Allocator& nodeAlloc, Allocator& leafAlloc, const size_t mode, const size_t threadID) 
    {
      /* stop toplevel recursion at some number of items */
      if (mode == CREATE_TOP_LEVEL && (current.size() <= topLevelItemThreshold || g_state->numBuildRecords >= MAX_TOP_LEVEL_BINS)) 
      {
        assert(g_state->numBuildRecords < NUM_TOP_LEVEL_BINS);
        g_state->buildRecords[g_state->numBuildRecords++] = current;
        return empty;
      }
      
      __aligned(64) SmallBuildRecord children[BVH4::N];
      
      /* create leaf node */
      if (unlikely(current.depth >= BVH4::maxBuildDepth || current.size() <= BVH4BuilderMorton::MORTON_LEAF_THRESHOLD)) {
        return createLeaf(current,nodeAlloc,leafAlloc,threadID);
      }
      
      /* fill all 4 children by always splitting the one with the largest surface area */
      size_t numChildren = 1;
      children[0] = current;
      
      do {
        
        /* find best child with largest bounding box area */
        int bestChild = -1;
        unsigned bestItems = 0;
        for (unsigned int i=0; i<numChildren; i++)
        {
          /* ignore leaves as they cannot get split */
          if (children[i].size() <= BVH4BuilderMorton::MORTON_LEAF_THRESHOLD)
            continue;
          
          /* remember child with largest area */
          if (children[i].size() > bestItems) { 
            bestItems = children[i].size();
            bestChild = i;
          }
        }
        if (bestChild == -1) break;
        
        /*! split best child into left and right child */
        __aligned(64) SmallBuildRecord left, right;
        split(children[bestChild],left,right);
                
        /* add new children left and right */
        left.depth = right.depth = current.depth+1;
        children[bestChild] = children[numChildren-1];
        children[numChildren-1] = left;
        children[numChildren+0] = right;
        numChildren++;
        
      } while (numChildren < BVH4::N);
      
      /* create leaf node if no split is possible */
      if (unlikely(numChildren == 1)) {
        BBox3f bounds; createSmallLeaf(this,current,leafAlloc,threadID,bounds); return bounds;
      }
      
      /* allocate node */
      Node* node = (Node*) nodeAlloc.malloc(sizeof(Node)); node->clear();
      *current.parent = bvh->encodeNode(node);
      
      /* recurse into each child */
      BBox3f bounds0 = empty;
      for (size_t i=0; i<numChildren; i++) 
      {
        children[i].parent = &node->child(i);
        
        if (children[i].size() <= BVH4BuilderMorton::MORTON_LEAF_THRESHOLD) {
          const BBox3f bounds = createLeaf(children[i],nodeAlloc,leafAlloc,threadID);
          bounds0.extend(bounds);
          node->set(i,bounds);
        } else {
          const BBox3f bounds = recurse(children[i],nodeAlloc,leafAlloc,mode,threadID);
          bounds0.extend(bounds);
          node->set(i,bounds);
        }
      }
      return bounds0;
    }
    
    __forceinline BBox3f BVH4BuilderMorton::leafBoundsTriangle1(NodeRef& ref)
    {
      BBox3f bounds = empty;
      size_t num; Triangle1* tri = (Triangle1*) ref.leaf(num);
      for (size_t i=0; i<num; i++) 
        bounds.extend(tri[i].bounds());
      return bounds;
    }
    
    __forceinline BBox3f BVH4BuilderMorton::leafBoundsTriangle4(NodeRef& ref)
    {
      BBox3f bounds = empty;
      size_t num; Triangle4* tri = (Triangle4*) ref.leaf(num);
      for (size_t i=0; i<num; i++) 
        bounds.extend(tri[i].bounds());
      return bounds;
    }
    
    __forceinline BBox3f BVH4BuilderMorton::leafBoundsTriangle1v(NodeRef& ref)
    {
      BBox3f bounds = empty;
      size_t num; Triangle1v* tri = (Triangle1v*) ref.leaf(num);
      for (size_t i=0; i<num; i++) 
        bounds.extend(tri[i].bounds());
      return bounds;
    }
    
    __forceinline BBox3f BVH4BuilderMorton::leafBoundsTriangle4v(NodeRef& ref)
    {
      BBox3f bounds = empty;
      size_t num; Triangle4v* tri = (Triangle4v*) ref.leaf(num);
      for (size_t i=0; i<num; i++) 
        bounds.extend(tri[i].bounds());
      return bounds;
    }
    
    __forceinline BBox3f BVH4BuilderMorton::node_bounds(NodeRef& ref) const
    {
      if (ref.isNode())
        return ref.node()->bounds();
      else
        return leafBounds(ref);
    }
    
    BBox3f BVH4BuilderMorton::refit_toplevel(NodeRef& ref) const
    { 
      /* stop here if we encounter a barrier */
      if (unlikely(ref.isBarrier())) {
        ref.clearBarrier();
        return node_bounds(ref);
      }
      
      /* return point bound for empty nodes */
      if (unlikely(ref == BVH4::emptyNode))
        return BBox3f(empty);
      
      /* this is a leaf node */
      if (unlikely(ref.isLeaf()))
	    return leafBounds(ref);
      
      /* recurse if this is an internal node */
      Node* node = ref.node();
      const BBox3f bounds0 = refit_toplevel(node->child(0));
      const BBox3f bounds1 = refit_toplevel(node->child(1));
      const BBox3f bounds2 = refit_toplevel(node->child(2));
      const BBox3f bounds3 = refit_toplevel(node->child(3));
      
      /* AOS to SOA transform */
      BBox<sse3f> bounds;
      transpose((ssef&)bounds0.lower,(ssef&)bounds1.lower,(ssef&)bounds2.lower,(ssef&)bounds3.lower,bounds.lower.x,bounds.lower.y,bounds.lower.z);
      transpose((ssef&)bounds0.upper,(ssef&)bounds1.upper,(ssef&)bounds2.upper,(ssef&)bounds3.upper,bounds.upper.x,bounds.upper.y,bounds.upper.z);
      
      /* set new bounds */
      node->lower_x = bounds.lower.x;
      node->lower_y = bounds.lower.y;
      node->lower_z = bounds.lower.z;
      node->upper_x = bounds.upper.x;
      node->upper_y = bounds.upper.y;
      node->upper_z = bounds.upper.z;
      
      /* return merged bounds */
      const float lower_x = reduce_min(bounds.lower.x);
      const float lower_y = reduce_min(bounds.lower.y);
      const float lower_z = reduce_min(bounds.lower.z);
      const float upper_x = reduce_max(bounds.upper.x);
      const float upper_y = reduce_max(bounds.upper.y);
      const float upper_z = reduce_max(bounds.upper.z);
      return BBox3f(Vec3fa(lower_x,lower_y,lower_z),
                    Vec3fa(upper_x,upper_y,upper_z));
    }

    void BVH4BuilderMorton::build_sequential_morton(size_t threadIndex, size_t threadCount) 
    {
      /* start measurement */
      double t0 = 0.0f;
      if (g_verbose >= 2) t0 = getSeconds();

      /* compute scene bounds */
      global_bounds = computeBounds();
      bvh->bounds = global_bounds.geometry;

      /* compute morton codes */
      if (mesh)
        computeMortonCodes(0,numPrimitives,mesh->id,0,morton);
      else
        computeMortonCodes(0,numPrimitives,0,0,morton);

      /* sort morton codes */
      std::sort(&morton[0],&morton[numPrimitives]); // FIMXE: use radix sort
      
#if defined(DEBUG)
      for (size_t i=1; i<numPrimitives; i++)
        assert(morton[i-1].code <= morton[i].code);
#endif	    
      
      SmallBuildRecord br;
      br.init(0,numPrimitives);
      br.parent = &bvh->root;
      br.depth = 1;
      
      /* perform first splits in single threaded mode */
      nodeAllocator.reset();
      primAllocator.reset();
      __aligned(64) Allocator nodeAlloc(nodeAllocator);
      __aligned(64) Allocator leafAlloc(primAllocator);
      recurse(br,nodeAlloc,leafAlloc,RECURSE,threadIndex);	    
            
      /* stop measurement */
      if (g_verbose >= 2) dt = getSeconds()-t0;
    }
    
    void BVH4BuilderMorton::build_parallel_morton(size_t threadIndex, size_t threadCount, size_t, size_t, TaskScheduler::Event* event) 
    {
      /* wait for all threads to enter */
      g_state->barrier.wait(threadIndex,threadCount);

      /* start measurement */
      double t0 = 0.0f;
      if (g_verbose >= 2) t0 = getSeconds();
      
      /* initialize thread state */
      initThreadState(threadIndex,threadCount);
      
	  /* let all thread except for control thread wait for work */
      if (threadIndex != 0) {
        scheduler.dispatchTaskMainLoop(threadIndex,threadCount);
        return;
      }
      
      /* compute scene bounds */
      global_bounds.reset();
      scheduler.dispatchTask( task_computeBounds, this, threadIndex, threadCount );
      bvh->bounds = global_bounds.geometry;

	  /* compute morton codes */
      scheduler.dispatchTask( task_computeMortonCodes, this, threadIndex, threadCount );   
      
      /* padding */
      MortonID32Bit* __restrict__ const dest = (MortonID32Bit*) nodeAllocator.data;
      for (size_t i=numPrimitives; i<( (numPrimitives+7)&(-8) ); i++) {
        dest[i].code  = 0xffffffff; 
        dest[i].index = 0;
      }
      
      /* sort morton codes */
      scheduler.dispatchTask( task_radixsort, this, threadIndex, threadCount );

#if defined(DEBUG)
      for (size_t i=1; i<numPrimitives; i++)
        assert(morton[i-1].code <= morton[i].code);
#endif	    
      
      /* build and extract top-level tree */
      g_state->numBuildRecords = 0;
      topLevelItemThreshold = (numPrimitives + threadCount-1)/(2*threadCount);
      
      SmallBuildRecord br;
      br.init(0,numPrimitives);
      br.parent = &bvh->root;
      br.depth = 1;
      
      /* perform first splits in single threaded mode */
      nodeAllocator.reset();
      primAllocator.reset();
      __aligned(64) Allocator nodeAlloc(nodeAllocator);
      __aligned(64) Allocator leafAlloc(primAllocator);
      recurse(br,nodeAlloc,leafAlloc,CREATE_TOP_LEVEL,threadIndex);	    
      
      /* sort all subtasks by size */
      insertionsort_decending<SmallBuildRecord>(g_state->buildRecords,g_state->numBuildRecords);
      
      /* build sub-trees */
      g_state->workStack.reset();
      scheduler.dispatchTask( task_recurseSubMortonTrees, this, threadIndex, threadCount );
      
      /* refit toplevel part of tree */
      refit_toplevel(bvh->root);
      
      /* end task */
      scheduler.releaseThreads(threadCount);
      
      /* stop measurement */
      if (g_verbose >= 2) dt = getSeconds()-t0;
    }
    
    Builder* BVH4BuilderMortonFast (void* bvh, BuildSource* source, Scene* scene, const size_t minLeafSize, const size_t maxLeafSize) {
      return new BVH4BuilderMorton((BVH4*)bvh,source,scene,NULL,minLeafSize,maxLeafSize);
    }
    
    Builder* BVH4BuilderMortonTriangleMeshFast (void* bvh, TriangleMeshScene::TriangleMesh* mesh, const size_t minLeafSize, const size_t maxLeafSize) {
      return new BVH4BuilderMorton((BVH4*)bvh,NULL,mesh->parent,mesh,minLeafSize,maxLeafSize);
    }
  }
}

  
  
