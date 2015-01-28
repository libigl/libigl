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
#include "bvh4i_builder_morton.h"
#include "bvh4i_statistics.h"
#include "bvh4i_builder_util.h"

#define BVH_NODE_PREALLOC_FACTOR 1.1f

//#define PROFILE

#define DBG(x)

namespace embree 
{
  namespace isa
  {
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    
    static double dt = 0.0f;
    
    BVH4iBuilderMorton::BVH4iBuilderMorton (BVH4i* bvh, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize)
    : bvh(bvh), source(source), scene((Scene*)geometry), topLevelItemThreshold(0), encodeShift(0), encodeMask(0), numBuildRecords(0), 
      morton(NULL), node(NULL), accel(NULL), numGroups(0), numPrimitives(0), numNodes(0), numAllocatedNodes(0)
    {
    }
    
    void BVH4iBuilderMorton::build(size_t threadIndex, size_t threadCount) 
    {
      if (g_verbose >= 2)
        std::cout << "building BVH4i with Morton builder ... " << std::flush;
      
      /* do some global inits first */
      init();

      LockStepTaskScheduler::init(TaskScheduler::getNumThreads());

#if defined(PROFILE)
      PING;
      DBG_PRINT(source->size());
      double dt_min = pos_inf;
      double dt_avg = 0.0f;
      double dt_max = neg_inf;
      for (size_t i=0; i<20; i++) 
      {
        TaskScheduler::executeTask(threadIndex,threadCount,_build_parallel_morton,this,TaskScheduler::getNumThreads(),"build_parallel_morton");
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

      TaskScheduler::executeTask(threadIndex,threadCount,_build_parallel_morton,this,TaskScheduler::getNumThreads(),"build_parallel_morton");
      
      if (g_verbose >= 2) {
        double perf = source->size()/dt*1E-6;
        std::cout << "[DONE] " << 1000.0f*dt << "ms (" << perf << " Mtris/s)" << std::endl;
        std::cout << BVH4iStatistics(bvh).str();
      }
#endif
    }
    
    void BVH4iBuilderMorton::init()
    {
      bvh->init();
      
      /* calculate total number of primrefs */
      size_t numPrimitivesOld = numPrimitives;
      numGroups     = source->groups();
      numPrimitives = source->size();
      
      size_t maxPrimsPerGroup = 0;
      for (size_t group=0; group<numGroups; group++) 
	{
	  if (unlikely(scene->get(group) == NULL)) continue;
	  if (scene->get(group)->type != TRIANGLE_MESH) continue;
	  const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(group);
	  if (unlikely(!mesh->isEnabled())) continue;

	  maxPrimsPerGroup = max(maxPrimsPerGroup,mesh->numTriangles);
	}
      
      /* calculate groupID, primID encoding */
      encodeShift = __bsr(maxPrimsPerGroup) + 1;
      assert( ((unsigned int)1 << encodeShift) > maxPrimsPerGroup);
      encodeMask = ((size_t)1 << encodeShift)-1;
      size_t maxGroups = ((size_t)1 << (31-encodeShift))-1;

      DBG(DBG_PRINT(numGroups));
      DBG(DBG_PRINT(numPrimitives));
      DBG(DBG_PRINT(encodeMask));
      DBG(DBG_PRINT(encodeShift));
      DBG(DBG_PRINT(maxGroups));

      if (maxPrimsPerGroup > encodeMask || numGroups > maxGroups)
      {
        DBG_PRINT(numGroups);
        DBG_PRINT(numPrimitives);
        DBG_PRINT(maxPrimsPerGroup);
        DBG_PRINT(encodeMask);
        DBG_PRINT(maxGroups);
        FATAL("ENCODING ERROR");
      }
      
      /* preallocate arrays */
      const size_t additional_size = 16 * CACHELINE_SIZE;
      
      if (numPrimitivesOld != numPrimitives)
      {
        /* free previously allocated memory */
        const size_t old_size_morton = numPrimitivesOld * sizeof(MortonID32Bit) + additional_size;
        const size_t old_size_node  = numAllocatedNodes * sizeof(BVHNode) + additional_size;
        const size_t old_size_accel = numPrimitivesOld * sizeof(Triangle1) + additional_size;
        if (morton) os_free(morton,old_size_morton);
        if (node  ) os_free(node  ,old_size_node);
        if (accel ) os_free(accel ,old_size_accel);
        
        /* allocated memory for primrefs,nodes, and accel */
        numAllocatedNodes = numPrimitives * BVH_NODE_PREALLOC_FACTOR;
        const size_t size_morton = numPrimitives * sizeof(MortonID32Bit) + additional_size;
        const size_t size_node   = numAllocatedNodes * sizeof(BVHNode) + additional_size;
        const size_t size_accel  = numPrimitives * sizeof(Triangle1) + additional_size;
        
        morton = (MortonID32Bit* ) os_malloc(size_morton); memset(morton,0,size_morton);
        node   = (BVHNode*)        os_malloc(size_node  ); memset(node  ,0,size_node);
        accel  = (Triangle1*)  os_malloc(size_accel ); memset(accel ,0,size_accel);	
      }
    }
    
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    
    void BVH4iBuilderMorton::initThreadState(const size_t threadID, const size_t numThreads)
    {
      const size_t startID = (threadID+0)*numPrimitives/numThreads;
      const size_t endID   = (threadID+1)*numPrimitives/numThreads;
      
      /* find first group containing startID */
      size_t group = 0, skipped = 0;
      for (; group<numGroups; group++) 
      {       
        Geometry* geom = scene->get(group);
        if (geom == NULL || geom->type != TRIANGLE_MESH) continue;
        TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
	if (unlikely(!mesh->isEnabled())) continue;

        const size_t numTriangles = mesh->numTriangles;	
        if (skipped + numTriangles > startID) break;
        skipped += numTriangles;
      }
      
      /* store start group and offset */
      thread_startGroup[threadID] = group;
      thread_startGroupOffset[threadID] = startID - skipped;
    }
    
    void BVH4iBuilderMorton::computeBounds(const size_t threadID, const size_t numThreads)
    {
      const size_t startID = (threadID+0)*numPrimitives/numThreads;
      const size_t endID   = (threadID+1)*numPrimitives/numThreads;
      
      Centroid_Scene_AABB bounds;
      bounds.reset();
      
      size_t currentID = startID;
      size_t offset = thread_startGroupOffset[threadID];
      
      for (size_t group = thread_startGroup[threadID]; group<numGroups; group++) 
      {       
        Geometry* geom = scene->get(group);
        if (geom == NULL || geom->type != TRIANGLE_MESH) continue;
        TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
	if (unlikely(!mesh->isEnabled())) continue;

        for (size_t i=offset; i<mesh->numTriangles && currentID < endID; i++, currentID++)	 
          bounds.extend(mesh->bounds(i));
        
        offset = 0;
        if (currentID == endID) break;
      }
      
      global_bounds.extend_atomic(bounds);    
    }
    
    void BVH4iBuilderMorton::computeMortonCodes(const size_t threadID, const size_t numThreads)
    {
      const size_t startID = (threadID+0)*numPrimitives/numThreads;
      const size_t endID   = (threadID+1)*numPrimitives/numThreads;
      
      /* store the morton codes temporarily in 'node' memory */
      MortonID32Bit* __restrict__ const dest = (MortonID32Bit*)node; 
      
      /* compute mapping from world space into 3D grid */
      const ssef base     = (ssef)global_bounds.centroid2.lower;
      const ssef diagonal = (ssef)global_bounds.centroid2.upper - (ssef)global_bounds.centroid2.lower;
      const ssef scale    = select(diagonal != 0, rcp(diagonal) * ssef(LATTICE_SIZE_PER_DIM * 0.99f),ssef(0.0f));
      
      size_t currentID = startID;
      size_t offset = thread_startGroupOffset[threadID];
      
      /* use SSE to calculate morton codes */
      size_t slots = 0;
      ssei ax = 0, ay = 0, az = 0, ai = 0;
      
      for (size_t group = thread_startGroup[threadID]; group<numGroups; group++) 
      {       
        Geometry* geom = scene->get(group);
        if (geom == NULL || geom->type != TRIANGLE_MESH) continue;
        TriangleMeshScene::TriangleMesh* mesh = (TriangleMeshScene::TriangleMesh*) geom;
	if (unlikely(!mesh->isEnabled())) continue;

        const size_t numTriangles = min(mesh->numTriangles-offset,endID-currentID);
        
        for (size_t i=0; i<numTriangles; i++)	  
        {
          const BBox3f b = mesh->bounds(offset+i);
          const ssef lower = (ssef)b.lower;
          const ssef upper = (ssef)b.upper;
          const ssef centroid2 = lower+upper;
          const ssei binID = ssei((centroid2-base)*scale);
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
            store4i(&dest[currentID-4],unpacklo(code,ai));
            store4i(&dest[currentID-2],unpackhi(code,ai));
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
    
    void BVH4iBuilderMorton::recreateMortonCodes(SmallBuildRecord& current) const
    {
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
        const ssef centroid2 = lower+upper;
        const ssei binID = ssei((centroid2-base)*scale);
        const unsigned int bx = extract<0>(binID);
        const unsigned int by = extract<1>(binID);
        const unsigned int bz = extract<2>(binID);
        const unsigned int code = bitInterleave(bx,by,bz); //FIXME use SSE/AVX
        morton[i].code  = code;
      }
      
      quicksort_insertionsort_ascending<MortonID32Bit,512>(morton,current.begin,current.end-1); 
      
#if defined(DEBUG)
      for (size_t i=current.begin; i<current.end-1; i++)
        assert(morton[i].code <= morton[i+1].code);
#endif	    
      
    }
    
    void BVH4iBuilderMorton::radixsort(const size_t threadID, const size_t numThreads)
    {
      const size_t startID = (threadID+0)*numPrimitives/numThreads;
      const size_t endID   = (threadID+1)*numPrimitives/numThreads;
      
      MortonID32Bit* __restrict__ mortonID[2];
      mortonID[0] = (MortonID32Bit*) morton; 
      mortonID[1] = (MortonID32Bit*) node;
      
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
        LockStepTaskScheduler::syncThreads(threadID,numThreads);
        
        /* calculate total number of items for each bucket */
        /* 32bit uint allows for better vectorization      */
        
        
#if 0
        __aligned(64) unsigned int total[RADIX_BUCKETS];
        for (size_t i=0; i<RADIX_BUCKETS; i++)
          total[i] = 0;
        
        for (size_t i=0; i<numThreads; i++)
          for (size_t j=0; j<RADIX_BUCKETS; j++)
            total[j] += radixCount[i][j];
#else
        __aligned(64) unsigned int inner_offset[RADIX_BUCKETS];
        
#define CHUNK 64
        
#pragma unroll(CHUNK)
        for (size_t i=0; i<RADIX_BUCKETS; i++)
          inner_offset[i] = 0;
        
        for (size_t j=0; j<RADIX_BUCKETS; j+=CHUNK)
          for (size_t i=0; i<threadID; i++)
#pragma unroll(CHUNK)
            for (size_t k=0;k<CHUNK;k++)
              inner_offset[j+k] += radixCount[i][j+k];
        
        __aligned(64) unsigned int total[RADIX_BUCKETS];
        
#pragma unroll(CHUNK)      
        for (size_t i=0; i<RADIX_BUCKETS; i++)
          total[i] = inner_offset[i];
        
        for (size_t j=0; j<RADIX_BUCKETS; j+=CHUNK)
          for (size_t i=threadID; i<numThreads; i++)
#pragma unroll(CHUNK)
            for (size_t k=0;k<CHUNK;k++)
              total[j+k] += radixCount[i][j+k];
        
#endif
        
        /* calculate start offset of each bucket */
        __aligned(64) unsigned int offset[RADIX_BUCKETS];
        offset[0] = 0;
        for (size_t i=1; i<RADIX_BUCKETS; i++)    
          offset[i] = offset[i-1] + total[i-1];
        
        /* calculate start offset of each bucket for this thread */
#if 0
        for (size_t j=0; j<RADIX_BUCKETS; j++)
          for (size_t i=0; i<threadID; i++)
            offset[j] += radixCount[i][j];
#else
        
#pragma unroll(32) 
        for (size_t j=0; j<RADIX_BUCKETS; j++)
          offset[j] += inner_offset[j];
#endif
        
        /* copy items into their buckets */
        for (size_t i=startID; i<endID; i++) {
          const unsigned int index = src[i].get(shift, mask);
          dst[offset[index]++] = src[i];
        }
        if (b < 2) LockStepTaskScheduler::syncThreads(threadID,numThreads);
      }
    }
    
    void BVH4iBuilderMorton::recurseSubMortonTrees(const size_t threadID, const size_t numThreads)
    {
      NodeAllocator alloc(atomicID,numAllocatedNodes);
      while (true)
      {
        const unsigned int taskID = LockStepTaskScheduler::taskCounter.inc();
        if (taskID >= numBuildRecords) break;
        
        recurse(buildRecords[taskID],alloc,RECURSE,threadID);
        
        /* mark toplevel of tree */
        node[buildRecords[taskID].parentID].upper.a = -1;
      }    
    }
    
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
      
      store4f(qptr->min_x,min_x);
      store4f(qptr->min_y,min_y);
      store4f(qptr->min_z,min_z);
      store4i(qptr->min_d,min_d);
      
      store4f(qptr->max_x,max_x);
      store4f(qptr->max_y,max_y);
      store4f(qptr->max_z,max_z);
      store4i(qptr->max_d,cast(max_a));       
    }
    
    void BVH4iBuilderMorton::convertToSOALayout(const size_t threadID, const size_t numThreads)
    {
      const size_t startID = (threadID+0)*numNodes/numThreads;
      const size_t endID   = (threadID+1)*numNodes/numThreads;
      
      BVHNode  * __restrict__  bptr = ( BVHNode*)node + startID*4;
      
      for (size_t n=startID; n<endID; n++, bptr+=4)
        convertToSOALayoutBlock4(bptr);
    }
    
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    __forceinline BBox3f BVH4iBuilderMorton::createSmallLeaf(SmallBuildRecord& current, NodeAllocator& alloc) const
    {
      ssef lower(pos_inf);
      ssef upper(neg_inf);
      size_t items = current.size();
      size_t start = current.begin;
      assert(items<=4);
      
      for (size_t i=0; i<items; i++) 
      {	
        const size_t index = morton[start+i].index;
        const size_t primID = index & encodeMask; 
        const size_t geomID = index >> encodeShift; 
        const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
        const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(primID);
        
        const ssef v0 = select(0x7,(ssef)mesh->vertex(tri.v[0]),zero);
        const ssef v1 = select(0x7,(ssef)mesh->vertex(tri.v[1]),zero);
        const ssef v2 = select(0x7,(ssef)mesh->vertex(tri.v[2]),zero);
        
        lower = min(lower,v0,v1,v2);
        upper = max(upper,v0,v1,v2);
        
        const ssef e1 = v0 - v1;
        const ssef e2 = v2 - v0;	     
        const ssef normal = cross(e1,e2);
        
        store4f_nt(&accel[start+i].v0,cast(insert<3>(cast(v0),primID)));
        store4f_nt(&accel[start+i].v1,cast(insert<3>(cast(v1),geomID)));
        store4f_nt(&accel[start+i].v2,cast(insert<3>(cast(v2),0)));
        store4f_nt(&accel[start+i].Ng,cast(insert<3>(cast(normal),0)));
      }
      node[current.parentID].lower = Vec3fa(lower);
      node[current.parentID].upper = Vec3fa(upper);
      node[current.parentID].createLeaf(start,items,items);
      BBox3f bounds;
      bounds.lower = Vec3fa(lower);
      bounds.upper = Vec3fa(upper);
      return bounds;
    }
    
    void BVH4iBuilderMorton::split_fallback(SmallBuildRecord& current, SmallBuildRecord& leftChild, SmallBuildRecord& rightChild) const
    {
      const unsigned int center = (current.begin + current.end)/2;
      leftChild.init(current.begin,center);
      rightChild.init(center,current.end);
    }
    
    BBox3f BVH4iBuilderMorton::createLeaf(SmallBuildRecord& current, NodeAllocator& alloc)
    {
#if defined(DEBUG)
      if (current.depth > BVH4i::maxBuildDepthLeaf) 
        throw std::runtime_error("ERROR: depth limit reached");
#endif
      
      /* create leaf for few primitives */
      if (current.size() <= MORTON_LEAF_THRESHOLD) {     
        return createSmallLeaf(current,alloc);
      }
      
      /* first split level */
      SmallBuildRecord record0, record1;
      split_fallback(current,record0,record1);
      
      /* second split level */
      SmallBuildRecord children[4];
      split_fallback(record0,children[0],children[1]);
      split_fallback(record1,children[2],children[3]);
      
      /* allocate next four nodes */
      size_t numChildren = 4;
      //const unsigned int currentIndex = allocNode(BVH4i::N);
      const size_t currentIndex = alloc.get(BVH4i::N);
      
      BBox3f bounds; 
      bounds = empty;
      /* recurse into each child */
      for (size_t i=0; i<numChildren; i++) {
        children[i].parentID = currentIndex+i;
        children[i].depth = current.depth+1;
        bounds.extend( createLeaf(children[i],alloc) );
      }
      
      node[current.parentID].lower = bounds.lower;
      node[current.parentID].upper = bounds.upper;
      node[current.parentID].createNode(currentIndex,numChildren);
      
      return bounds;
    }  
    
    __forceinline bool BVH4iBuilderMorton::split(SmallBuildRecord& current,
                                                 SmallBuildRecord& left,
                                                 SmallBuildRecord& right) const
    {
      /* mark as leaf if leaf threshold reached */
      if (unlikely(current.size() <= BVH4iBuilderMorton::MORTON_LEAF_THRESHOLD)) {
        //current.createLeaf();
        return false;
      }
      
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
          return true;
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
      return true;
    }
    
    BBox3f BVH4iBuilderMorton::recurse(SmallBuildRecord& current, 
                                       NodeAllocator& alloc,
                                       const size_t mode, 
                                       const size_t threadID) 
    {
      /* stop toplevel recursion at some number of items */
      if (mode == CREATE_TOP_LEVEL && current.size() <= topLevelItemThreshold) {
        buildRecords[numBuildRecords++] = current; // FIXME: can overflow
        return empty;
      }
      
      __aligned(64) SmallBuildRecord children[BVH4i::N];
      
      /* create leaf node */
      if (unlikely(current.size() <= BVH4iBuilderMorton::MORTON_LEAF_THRESHOLD)) {
        return createSmallLeaf(current,alloc);
      }
      if (unlikely(current.depth >= BVH4i::maxBuildDepth)) {
        return createLeaf(current,alloc); 
      }
      
      /* fill all 4 children by always splitting the one with the largest surface area */
      size_t numChildren = 1;
      children[0] = current;
      
      do {
        
        /* find best child with largest number of items*/
        int bestChild = -1;
        unsigned bestItems = 0;
        for (unsigned int i=0; i<numChildren; i++)
        {
          /* ignore leaves as they cannot get split */
          if (children[i].size() <= BVH4iBuilderMorton::MORTON_LEAF_THRESHOLD)
            continue;
          
          /* remember child with largest number of items */
          if (children[i].size() > bestItems) { 
            bestItems = children[i].size();
            bestChild = i;
          }
        }
        if (bestChild == -1) break;
        
        /*! split best child into left and right child */
        __aligned(64) SmallBuildRecord left, right;
        if (!split(children[bestChild],left,right))
          continue;
        
        /* add new children left and right */
        left.depth = right.depth = current.depth+1;
        children[bestChild] = children[numChildren-1];
        children[numChildren-1] = left;
        children[numChildren+0] = right;
        numChildren++;
        
      } while (numChildren < BVH4i::N);
      
      /* create leaf node if no split is possible */
      if (unlikely(numChildren == 1)) {
        return createSmallLeaf(current,alloc);
      }
      
      /* allocate next four nodes */
      const size_t currentIndex = alloc.get(BVH4i::N);
      
      /* recurse into each child */
      BBox3f bounds;
      bounds = empty;
      for (size_t i=0; i<numChildren; i++) 
      {
        children[i].parentID = currentIndex+i;
        
        if (children[i].size() <= BVH4iBuilderMorton::MORTON_LEAF_THRESHOLD)
          bounds.extend( createLeaf(children[i],alloc) );
        else
          bounds.extend( recurse(children[i],alloc,mode,threadID) );
      }
      node[current.parentID].lower = bounds.lower;
      node[current.parentID].upper = bounds.upper;
      node[current.parentID].createNode(currentIndex,numChildren);
      
      /* init unused nodes */
      const avxf init_node = load8f((float*)initQBVHNode);
      for (size_t i=numChildren; i<BVH4i::N; i++)
        store8f_nt((float*)&node[currentIndex+i],init_node);
      
      return bounds;
    }
    
    
    void BVH4iBuilderMorton::refit(const size_t index) const
    {    
      BVHNode& entry = node[index];
      
      if (unlikely(entry.isLeaf()))
        return;
      
      const size_t children = entry.firstChildID();
      const size_t items    = entry.items();
      BVHNode* next = &node[children+0];
      
      ssef lower(pos_inf);
      ssef upper(neg_inf);
      const int e0 = entry.lower.a;
      const int e1 = entry.upper.a;
      
      for (size_t i=0; i<items; i++) 
      {
        const size_t childIndex = children + i;	    	    
        if (!next[i].isLeaf())
          refit(childIndex);
        
        lower = min(lower,(ssef)next[i].lower);
        upper = max(upper,(ssef)next[i].upper);
      }      
      
      entry.lower = Vec3fa(lower);
      entry.upper = Vec3fa(upper);
      entry.lower.a = e0;
      entry.upper.a = e1;
    }    
    
    void BVH4iBuilderMorton::refit_toplevel(const size_t index) const
    {    
      BVHNode& entry = node[index];
      
      if (entry.upper.a == -1 || entry.isLeaf())    
        return;
      
      const unsigned int children = entry.firstChildID();
      BVHNode* next = &node[children+0];
      
      const unsigned int items = entry.items();
      
      ssef lower(pos_inf);
      ssef upper(neg_inf);
      const int e0 = entry.lower.a;
      const int e1 = entry.upper.a;
      
      for (unsigned int i=0; i<items; i++) 
      {
        const unsigned int childIndex = children + i;	    	    
        if (!next[i].isLeaf())
          refit_toplevel(childIndex);
        
        lower = min(lower,(ssef)next[i].lower);
        upper = max(upper,(ssef)next[i].upper);
      }      
      
      entry.lower = Vec3fa(lower);
      entry.upper = Vec3fa(upper);
      entry.lower.a = e0;
      entry.upper.a = e1;
    }
    
    void BVH4iBuilderMorton::build_main (const size_t threadIndex, const size_t threadCount)
    { 
      /* compute scene bounds */
      global_bounds.reset();
      LockStepTaskScheduler::dispatchTask( task_computeBounds, this, threadIndex, threadCount );
      
      /* compute morton codes */
      LockStepTaskScheduler::dispatchTask( task_computeMortonCodes, this, threadIndex, threadCount );   
      
      /* padding */
      MortonID32Bit* __restrict__ const dest = (MortonID32Bit*)node;
      for (size_t i=numPrimitives; i<( (numPrimitives+7)&(-8) ); i++) {
        dest[i].code  = 0xffffffff; 
        dest[i].index = 0;
      }
      
      /* sort morton codes */
      LockStepTaskScheduler::dispatchTask( task_radixsort, this, threadIndex, threadCount );
      
#if defined(DEBUG)
      for (size_t i=1; i<numPrimitives; i++)
        assert(morton[i-1].code <= morton[i].code);
#endif	    
      
      /* build and extract top-level tree */
      numBuildRecords = 0;
      atomicID.reset(BVH4i::N);
      node[0].lower = global_bounds.geometry.lower;
      node[0].upper = global_bounds.geometry.upper;
      topLevelItemThreshold = (numPrimitives + threadCount-1)/(2*threadCount);
      
      SmallBuildRecord br;
      br.init(0,numPrimitives);
      br.parentID = 0;
      br.depth = 1;
      
      /* perform first splits in single threaded mode */
      NodeAllocator alloc(atomicID,numAllocatedNodes);
      recurse(br,alloc,CREATE_TOP_LEVEL,threadIndex);	    
      
      /* sort all subtasks by size */
      insertionsort_decending<SmallBuildRecord>(buildRecords,numBuildRecords);
      
      /* build sub-trees */
      LockStepTaskScheduler::dispatchTask( task_recurseSubMortonTrees, this, threadIndex, threadCount );
      numNodes = atomicID >> 2;
      
      /* refit toplevel part of tree */
      refit_toplevel(0);
    }
    
    void BVH4iBuilderMorton::build_parallel_morton(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event) 
    {
      /* init lock step task scheduler */
      if (threadIndex == 0) 
        LockStepTaskScheduler::init(TaskScheduler::getNumThreads()); 
      
      /* start measurement */
      double t0 = 0.0f;
      if (g_verbose >= 2) t0 = getSeconds();
      
      /* initialize thread state */
      initThreadState(threadIndex,threadCount);
      
      /* let all thread except for control thread wait for work */
      if (threadIndex != 0) {
        LockStepTaskScheduler::dispatchTaskMainLoop(threadIndex,threadCount);
        return;
      }
      
      /* performs build of tree */
      build_main(threadIndex,taskCount);
      
      /* convert to optimized layout */
      bvh->accel = accel;
      bvh->qbvh  = node;
      LockStepTaskScheduler::dispatchTask( task_convertToSOALayout, this, threadIndex, threadCount );
      
      /* set root and bounding box */
      const QBVHNode* const qbvh  = (QBVHNode*)bvh->qbvh;
      bvh->root = qbvh[0].min_d[0]; 
      bvh->bounds = global_bounds.geometry;
      
      /* end task */
      LockStepTaskScheduler::releaseThreads(threadCount);
      
      /* stop measurement */
      if (g_verbose >= 2) dt = getSeconds()-t0;
    }
 
    Builder* BVH4iTriangle1BuilderMorton (void* bvh, BuildSource* source, Scene* scene, const size_t minLeafSize, const size_t maxLeafSize) {
      return new BVH4iBuilderMorton((BVH4i*)bvh,source,scene,minLeafSize,maxLeafSize);
    }
  }
}
  
  
