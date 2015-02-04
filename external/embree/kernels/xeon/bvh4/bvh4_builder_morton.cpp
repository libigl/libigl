// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
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
#include "geometry/triangle8.h"
#include "geometry/triangle1v.h"
#include "geometry/triangle4v.h"
#include "geometry/triangle4i.h"

#include "sys/tasklogger.h"

#define DBG(x) 

namespace embree 
{
  namespace isa
  {
    static double dt = 0.0f;

    BVH4BuilderMorton::BVH4BuilderMorton (BVH4* bvh, Scene* scene, TriangleMesh* mesh, size_t listMode, size_t logBlockSize, bool needVertices, size_t primBytes, const size_t minLeafSize, const size_t maxLeafSize)
      : bvh(bvh), state(nullptr), scheduler(&scene->lockstep_scheduler), scene(scene), mesh(mesh), listMode(listMode), logBlockSize(logBlockSize), needVertices(needVertices), primBytes(primBytes), minLeafSize(minLeafSize), maxLeafSize(maxLeafSize),
	topLevelItemThreshold(0), encodeShift(0), encodeMask(-1), morton(NULL), bytesMorton(0), numGroups(0), numPrimitives(0), numAllocatedPrimitives(0), numAllocatedNodes(0)
    {
      needAllThreads = true;
      if (mesh) needAllThreads = mesh->numTriangles > 50000;
    }
    
    BVH4Triangle1BuilderMorton::BVH4Triangle1BuilderMorton (BVH4* bvh, Scene* scene, size_t listMode)
      : BVH4BuilderMorton(bvh,scene,NULL,listMode,0,false,sizeof(Triangle1),4,inf) {}

    BVH4Triangle4BuilderMorton::BVH4Triangle4BuilderMorton (BVH4* bvh, Scene* scene, size_t listMode)
      : BVH4BuilderMorton(bvh,scene,NULL,listMode,2,false,sizeof(Triangle4),4,inf) {}

#if defined(__AVX__)
    BVH4Triangle8BuilderMorton::BVH4Triangle8BuilderMorton (BVH4* bvh, Scene* scene, size_t listMode)
      : BVH4BuilderMorton(bvh,scene,NULL,listMode,3,false,sizeof(Triangle8),8,inf) {}
#endif
    
    BVH4Triangle1vBuilderMorton::BVH4Triangle1vBuilderMorton (BVH4* bvh, Scene* scene, size_t listMode)
      : BVH4BuilderMorton(bvh,scene,NULL,listMode,0,false,sizeof(Triangle1v),4,inf) {}

    BVH4Triangle4vBuilderMorton::BVH4Triangle4vBuilderMorton (BVH4* bvh, Scene* scene, size_t listMode)
      : BVH4BuilderMorton(bvh,scene,NULL,listMode,2,false,sizeof(Triangle4v),4,inf) {}

    BVH4Triangle4iBuilderMorton::BVH4Triangle4iBuilderMorton (BVH4* bvh, Scene* scene, size_t listMode)
      : BVH4BuilderMorton(bvh,scene,NULL,listMode,2,true,sizeof(Triangle4i),4,inf) {}

    BVH4Triangle1BuilderMorton::BVH4Triangle1BuilderMorton (BVH4* bvh, TriangleMesh* mesh, size_t listMode)
      : BVH4BuilderMorton(bvh,mesh->parent,mesh,listMode,0,false,sizeof(Triangle1),4,inf) {}

    BVH4Triangle4BuilderMorton::BVH4Triangle4BuilderMorton (BVH4* bvh, TriangleMesh* mesh, size_t listMode)
      : BVH4BuilderMorton(bvh,mesh->parent,mesh,listMode,2,false,sizeof(Triangle4),4,inf) {}

#if defined(__AVX__)
    BVH4Triangle8BuilderMorton::BVH4Triangle8BuilderMorton (BVH4* bvh, TriangleMesh* mesh, size_t listMode)
      : BVH4BuilderMorton(bvh,mesh->parent,mesh,listMode,3,false,sizeof(Triangle8),8,inf) {}
#endif
    
    BVH4Triangle1vBuilderMorton::BVH4Triangle1vBuilderMorton (BVH4* bvh, TriangleMesh* mesh, size_t listMode)
      : BVH4BuilderMorton(bvh,mesh->parent,mesh,listMode,0,false,sizeof(Triangle1v),4,inf) {}

    BVH4Triangle4vBuilderMorton::BVH4Triangle4vBuilderMorton (BVH4* bvh, TriangleMesh* mesh, size_t listMode)
      : BVH4BuilderMorton(bvh,mesh->parent,mesh,listMode,2,false,sizeof(Triangle4v),4,inf) {}

    BVH4Triangle4iBuilderMorton::BVH4Triangle4iBuilderMorton (BVH4* bvh, TriangleMesh* mesh, size_t listMode)
      : BVH4BuilderMorton(bvh,mesh->parent,mesh,listMode,2,true,sizeof(Triangle4i),4,inf) {}
        
    BVH4BuilderMorton::~BVH4BuilderMorton () 
    {
      if (morton) os_free(morton,bytesMorton);
      bvh->alloc.shrink();
    }
    
    void BVH4BuilderMorton::build(size_t threadIndex, size_t threadCount) 
    {
      if (g_verbose >= 2)
        std::cout << "building BVH4<" << bvh->primTy.name << "> with " << TOSTRING(isa) << "::BVH4BuilderMorton ... " << std::flush;
      
      /* calculate size of scene */
      size_t numPrimitivesOld = numPrimitives;
      if (mesh) numPrimitives = mesh->numTriangles;
      else      numPrimitives = scene->numTriangles;
      
      bvh->numPrimitives = numPrimitives;
      numGroups = scene->size();
      size_t maxPrimsPerGroup = 0;
      if (mesh) 
        maxPrimsPerGroup = numPrimitives;
      else {
	for (size_t i=0; i<numGroups; i++) // FIXME: encoding unsafe
        {
          Geometry* geom = scene->get(i);
          if (!geom || geom->type != TRIANGLE_MESH) continue;
          TriangleMesh* mesh = (TriangleMesh*) geom;
	  if (mesh->numTimeSteps != 1) continue;
	  maxPrimsPerGroup = max(maxPrimsPerGroup,mesh->numTriangles);
        }

        /* calculate groupID, primID encoding */
        encodeShift = __bsr(maxPrimsPerGroup) + 1;
        encodeMask = ((size_t)1 << encodeShift)-1;
        size_t maxGroups = ((size_t)1 << (31-encodeShift))-1;
      
        if (maxPrimsPerGroup > encodeMask || numGroups > maxGroups) 
          THROW_RUNTIME_ERROR("encoding error in morton builder");
      }
      
      /* preallocate arrays */
      if (numPrimitivesOld != numPrimitives)
      {
	bvh->init(sizeof(BVH4::Node),numPrimitives,threadCount);
        if (morton) os_free(morton,bytesMorton);
	bytesMorton = ((numPrimitives+7)&(-8)) * sizeof(MortonID32Bit);
        morton = (MortonID32Bit* ) os_malloc(bytesMorton); memset(morton,0,bytesMorton);
      }
         
      if (needAllThreads) 
      {
        state.reset(new MortonBuilderState);
	//size_t numActiveThreads = threadCount;
	size_t numActiveThreads = min(threadCount,getNumberOfCores());
	build_parallel_morton(threadIndex,numActiveThreads,0,1);
        state.reset(NULL);
      } else {
        build_sequential_morton(threadIndex,threadCount);
      }

      if (g_verbose >= 2) {
        std::cout << "[DONE] " << 1000.0f*dt << "ms (" << numPrimitives/dt*1E-6 << " Mtris/s)" << std::endl;
        std::cout << BVH4Statistics(bvh).str();
      }

      /* benchmark mode */
      if (g_benchmark) {
	BVH4Statistics stat(bvh);
	std::cout << "BENCHMARK_BUILD " << dt << " " << double(numPrimitives)/dt << " " << stat.sah() << " " << stat.bytesUsed() << std::endl;
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
        state->startGroup[threadID] = mesh->id;
        state->startGroupOffset[threadID] = startID;
      }
      else
      {
        /* find first group containing startID */
        size_t group = 0, skipped = 0;
        for (; group<numGroups; group++) 
        {       
          Geometry* geom = scene->get(group);
          if (!geom || geom->type != TRIANGLE_MESH) continue;
          TriangleMesh* mesh = (TriangleMesh*) geom;
          if (mesh->numTimeSteps != 1) continue;
          const size_t numTriangles = mesh->numTriangles;
          if (skipped + numTriangles > startID) break; 
         skipped += numTriangles;
        }
        
        /* store start group and offset */
        state->startGroup[threadID] = group;
        state->startGroupOffset[threadID] = startID - skipped;
      }
    }

    CentGeomBBox3fa BVH4BuilderMorton::computeBounds()
    {
      CentGeomBBox3fa bounds; bounds.reset();

      if (mesh) 
      {
        for (size_t i=0; i<mesh->numTriangles; i++) {
	  const BBox3fa b = mesh->bounds(i);
	  if (!inFloatRange(b)) continue;
          bounds.extend(b);
	}
      }
      else
      {
        for (size_t group=0; group<numGroups; group++) 
        {       
          Geometry* geom = scene->get(group);
          if (!geom || geom->type != TRIANGLE_MESH) continue;
          TriangleMesh* mesh = (TriangleMesh*) geom;
          if (mesh->numTimeSteps != 1) continue;
          for (size_t i=0; i<mesh->numTriangles; i++) {
	    const BBox3fa b = mesh->bounds(i);
	    if (!inFloatRange(b)) continue;
	    bounds.extend(b);
	  }
        }
      }
      return bounds;
    }

    void BVH4BuilderMorton::computeBounds(const size_t threadID, const size_t numThreads)
    {
      initThreadState(threadID,numThreads);

      const ssize_t start = (threadID+0)*numPrimitives/numThreads;
      const ssize_t end   = (threadID+1)*numPrimitives/numThreads;
      
      CentGeomBBox3fa bounds; bounds.reset();
      
      if (mesh) 
      {
        for (size_t i=start; i<end; i++) {
	  const BBox3fa b = mesh->bounds(i);
	  if (!inFloatRange(b)) continue;
          bounds.extend(b);
	}
      }
      else
      {
	for (ssize_t cur=0, i=0; i<ssize_t(scene->size()); i++) 
	{
	  TriangleMesh* geom = (TriangleMesh*) scene->get(i);
          if (geom == NULL) continue;
	  if (geom->type != TRIANGLE_MESH || geom->numTimeSteps != 1 || !geom->isEnabled()) continue;
	  ssize_t gstart = 0;
	  ssize_t gend = geom->numTriangles;
	  ssize_t s = max(start-cur,gstart);
	  ssize_t e = min(end  -cur,gend  );
	  for (ssize_t j=s; j<e; j++) {
	    const BBox3fa b = geom->bounds(j);
	    if (!inFloatRange(b)) continue;
	    bounds.extend(b);
	  }
	  cur += geom->numTriangles;
	  if (cur >= end) break;  
        }
      }
      global_bounds.extend_atomic(bounds);    
    }

    void BVH4BuilderMorton::computeMortonCodes(const size_t startID, const size_t endID, size_t& destID,
                                               const size_t startGroup, const size_t startOffset, 
                                               MortonID32Bit* __restrict__ const dest)
    {
      /* compute mapping from world space into 3D grid */
      const ssef base  = (ssef)global_bounds.centBounds.lower;
      const ssef diag  = (ssef)global_bounds.centBounds.upper - (ssef)global_bounds.centBounds.lower;
      const ssef scale = select(diag > ssef(1E-19f), rcp(diag) * ssef(LATTICE_SIZE_PER_DIM * 0.99f),ssef(0.0f));
      
      size_t currentID = destID;
      size_t offset = startOffset;
      
      /* use SSE to calculate morton codes */
      size_t slots = 0;
      ssei ax = 0, ay = 0, az = 0, ai = 0;
      
      for (size_t group = startGroup; group<numGroups; group++) 
      {       
        Geometry* geom = scene->get(group);
        if (!geom || !geom->isEnabled() || geom->type != TRIANGLE_MESH) continue;
        TriangleMesh* mesh = (TriangleMesh*) geom;
        if (mesh->numTimeSteps != 1) continue;
        const size_t numTriangles = min(mesh->numTriangles-offset,endID-currentID);
        
        for (size_t i=0; i<numTriangles; i++)	  
        {
          const BBox3fa b = mesh->bounds(offset+i);
          const ssef lower = (ssef)b.lower;
          const ssef upper = (ssef)b.upper;
          const ssef centroid = lower+upper;
          const ssei binID = ssei((centroid-base)*scale);
          unsigned int index = offset+i;
          if (this->mesh == NULL) index |= group << encodeShift;
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
      destID = currentID - destID;
    }
    
    void BVH4BuilderMorton::computeMortonCodes(const size_t threadID, const size_t numThreads)
    {      
      const size_t startID = (threadID+0)*numPrimitives/numThreads;
      const size_t endID   = (threadID+1)*numPrimitives/numThreads;
      
      /* store the morton codes temporarily in 'node' memory */
      MortonID32Bit* __restrict__ const dest = (MortonID32Bit*)bvh->alloc.base();
      computeMortonCodes(startID,endID,state->dest[threadID],state->startGroup[threadID],state->startGroupOffset[threadID],dest);
    }
    
    void BVH4BuilderMorton::recreateMortonCodes(BuildRecord& current) const
    {
      assert(current.size() > 4);
      CentGeomBBox3fa global_bounds;
      global_bounds.reset();
      
      for (size_t i=current.begin; i<current.end; i++)
      {
        const size_t index  = morton[i].index;
        const size_t primID = index & encodeMask; 
        const size_t geomID = index >> encodeShift; 
        const TriangleMesh* mesh = this->mesh ? this->mesh : scene->getTriangleMesh(geomID);
        global_bounds.extend(mesh->bounds(primID));
      }
      
      /* compute mapping from world space into 3D grid */
      const ssef base  = (ssef)global_bounds.centBounds.lower;
      const ssef diag  = (ssef)global_bounds.centBounds.upper - (ssef)global_bounds.centBounds.lower;
      const ssef scale = select(diag > ssef(1E-19f), rcp(diag) * ssef(LATTICE_SIZE_PER_DIM * 0.99f),ssef(0.0f));
      
      for (size_t i=current.begin; i<current.end; i++)
      {
        const size_t index  = morton[i].index;
        const size_t primID = index & encodeMask; 
        const size_t geomID = index >> encodeShift; 
        const TriangleMesh* mesh = this->mesh ? this->mesh : scene->getTriangleMesh(geomID);
        const BBox3fa b = mesh->bounds(primID);
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
      std::sort(morton+current.begin,morton+current.end);
      
#if defined(DEBUG)
      for (size_t i=current.begin; i<current.end-1; i++)
        assert(morton[i].code <= morton[i+1].code);
#endif	    
    }
    
    void BVH4BuilderMorton::radixsort(const size_t threadID, const size_t numThreads)
    {
      const size_t startID = (threadID+0)*numPrimitives/numThreads;
      const size_t endID   = (threadID+1)*numPrimitives/numThreads;
      
      MortonID32Bit* __restrict__ mortonID[2];
      mortonID[0] = (MortonID32Bit*) morton; 
      mortonID[1] = (MortonID32Bit*) bvh->alloc.base();
      MortonBuilderState::ThreadRadixCountTy* radixCount = state->radixCount;
      
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
        //TaskScheduler::syncThreads(threadID,numThreads);
	barrier.wait(threadID,numThreads);
        
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
        if (b < 2) barrier.wait(threadID,numThreads);
	  //TaskScheduler::syncThreads(threadID,numThreads);
      }
    }
    
    void BVH4BuilderMorton::recurseSubMortonTrees(const size_t threadID, const size_t numThreads)
    {
      __aligned(64) Allocator nodeAlloc(&bvh->alloc);
      __aligned(64) Allocator leafAlloc(&bvh->alloc);
      while (true)
      {
        const unsigned int taskID = atomic_add(&state->taskCounter,1);
        if (taskID >= state->buildRecords.size()) break;
	recurse(state->buildRecords[taskID],nodeAlloc,leafAlloc,RECURSE,threadID);
        state->buildRecords[taskID].parent->setBarrier();
      }
      _mm_sfence(); // make written leaves globally visible
    }
    
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    
    void BVH4Triangle1BuilderMorton::createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3fa& box_o)
    {
      ssef lower(pos_inf);
      ssef upper(neg_inf);
      size_t items = current.size();
      size_t start = current.begin;
      
      /* allocate leaf node */
      Triangle1* accel = (Triangle1*) leafAlloc.malloc(items*sizeof(Triangle1));
      *current.parent = bvh->encodeLeaf((char*)accel,listMode ? listMode : items);
      
      for (size_t i=0; i<items; i++) 
      {	
        const size_t index = morton[start+i].index;
        const size_t primID = index & encodeMask; 
        const size_t geomID = this->mesh ? this->mesh->id : (index >> encodeShift); 
        const TriangleMesh* mesh = scene->getTriangleMesh(geomID);
        const TriangleMesh::Triangle& tri = mesh->triangle(primID);
        
        const ssef v0 = select(0x7,(ssef)mesh->vertex(tri.v[0]),zero);
        const ssef v1 = select(0x7,(ssef)mesh->vertex(tri.v[1]),zero);
        const ssef v2 = select(0x7,(ssef)mesh->vertex(tri.v[2]),zero);
        
        lower = min(lower,v0,v1,v2);
        upper = max(upper,v0,v1,v2);
        
        const ssef e1 = v0 - v1;
        const ssef e2 = v2 - v0;	     
        const ssef normal = cross(e1,e2);
	const bool last = listMode && (i==(items-1));
        
        store4f_nt(&accel[i].v0,cast(insert<3>(cast(v0),primID | (last << 31))));
        store4f_nt(&accel[i].v1,cast(insert<3>(cast(v1),geomID)));
        store4f_nt(&accel[i].v2,cast(insert<3>(cast(v2),mesh->mask)));
        store4f_nt(&accel[i].Ng,cast(insert<3>(cast(normal),0)));
      }
      box_o = BBox3fa((Vec3fa)lower,(Vec3fa)upper);
    }
    
    void BVH4Triangle4BuilderMorton::createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3fa& box_o)
    {
      ssef lower(pos_inf);
      ssef upper(neg_inf);
      size_t items = current.size();
      size_t start = current.begin;
      assert(items<=4);
      
      /* allocate leaf node */
      Triangle4* accel = (Triangle4*) leafAlloc.malloc(sizeof(Triangle4));
      *current.parent = bvh->encodeLeaf((char*)accel,listMode ? listMode : 1);
      
      ssei vgeomID = -1, vprimID = -1, vmask = -1;
      sse3f v0 = zero, v1 = zero, v2 = zero;
      
      for (size_t i=0; i<items; i++)
      {
        const size_t index = morton[start+i].index;
        const size_t primID = index & encodeMask; 
        const size_t geomID = this->mesh ? this->mesh->id : (index >> encodeShift); 
        const TriangleMesh* mesh = scene->getTriangleMesh(geomID);
        const TriangleMesh::Triangle& tri = mesh->triangle(primID);
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
      Triangle4::store_nt(accel,Triangle4(v0,v1,v2,vgeomID,vprimID,vmask,listMode));
      box_o = BBox3fa((Vec3fa)lower,(Vec3fa)upper);
    }

#if defined(__AVX__)
    void BVH4Triangle8BuilderMorton::createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3fa& box_o)
    {
      ssef lower(pos_inf);
      ssef upper(neg_inf);
      size_t items = current.size();
      size_t start = current.begin;
      assert(items<=8);
      
      /* allocate leaf node */
      Triangle8* accel = (Triangle8*) leafAlloc.malloc(sizeof(Triangle8));
      *current.parent = bvh->encodeLeaf((char*)accel,listMode ? listMode : 1);
      
      avxi vgeomID = -1, vprimID = -1, vmask = -1;
      avx3f v0 = zero, v1 = zero, v2 = zero;
      
      for (size_t i=0; i<items; i++)
      {
        const size_t index = morton[start+i].index;
        const size_t primID = index & encodeMask; 
        const size_t geomID = this->mesh ? this->mesh->id : (index >> encodeShift); 
        const TriangleMesh* mesh = scene->getTriangleMesh(geomID);
        const TriangleMesh::Triangle& tri = mesh->triangle(primID);
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
      new (accel) Triangle8(v0,v1,v2,vgeomID,vprimID,vmask,listMode); // FIXME: use storent
      box_o = BBox3fa((Vec3fa)lower,(Vec3fa)upper);
    }
#endif
    
    void BVH4Triangle1vBuilderMorton::createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3fa& box_o)
    {
      ssef lower(pos_inf);
      ssef upper(neg_inf);
      size_t items = current.size();
      size_t start = current.begin;
      
      /* allocate leaf node */
      Triangle1v* accel = (Triangle1v*) leafAlloc.malloc(items*sizeof(Triangle1v));
      *current.parent = bvh->encodeLeaf((char*)accel,listMode ? listMode : items);
      
      for (size_t i=0; i<items; i++) 
      {	
        const size_t index = morton[start+i].index;
        const size_t primID = index & encodeMask; 
        const size_t geomID = this->mesh ? this->mesh->id : (index >> encodeShift); 
        const TriangleMesh* mesh = scene->getTriangleMesh(geomID);
        const TriangleMesh::Triangle& tri = mesh->triangle(primID);
        
        const ssef v0 = select(0x7,(ssef)mesh->vertex(tri.v[0]),zero);
        const ssef v1 = select(0x7,(ssef)mesh->vertex(tri.v[1]),zero);
        const ssef v2 = select(0x7,(ssef)mesh->vertex(tri.v[2]),zero);
        
        lower = min(lower,v0,v1,v2);
        upper = max(upper,v0,v1,v2);
        
        const ssef e1 = v0 - v1;
        const ssef e2 = v2 - v0;	     
        const ssef normal = cross(e1,e2);
        const bool last = listMode && (i==(items-1));

        store4f_nt(&accel[i].v0,cast(insert<3>(cast(v0),primID | (last << 31))));
        store4f_nt(&accel[i].v1,cast(insert<3>(cast(v1),geomID)));
        store4f_nt(&accel[i].v2,cast(insert<3>(cast(v2),mesh->mask)));
      }
      box_o = BBox3fa((Vec3fa)lower,(Vec3fa)upper);
    }
    
    void BVH4Triangle4vBuilderMorton::createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3fa& box_o)
    {
      ssef lower(pos_inf);
      ssef upper(neg_inf);
      size_t items = current.size();
      size_t start = current.begin;
      assert(items<=4);
      
      /* allocate leaf node */
      Triangle4v* accel = (Triangle4v*) leafAlloc.malloc(sizeof(Triangle4v));
      *current.parent = bvh->encodeLeaf((char*)accel,listMode ? listMode : 1);
      
      ssei vgeomID = -1, vprimID = -1, vmask = -1;
      sse3f v0 = zero, v1 = zero, v2 = zero;
      
      for (size_t i=0; i<items; i++)
      {
        const size_t index = morton[start+i].index;
        const size_t primID = index & encodeMask; 
        const size_t geomID = this->mesh ? this->mesh->id : (index >> encodeShift); 
        const TriangleMesh* mesh = scene->getTriangleMesh(geomID);
        const TriangleMesh::Triangle& tri = mesh->triangle(primID);
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
      Triangle4v::store_nt(accel,Triangle4v(v0,v1,v2,vgeomID,vprimID,vmask,listMode));
      box_o = BBox3fa((Vec3fa)lower,(Vec3fa)upper);
    }

    void BVH4Triangle4iBuilderMorton::createSmallLeaf(BuildRecord& current, Allocator& leafAlloc, size_t threadID, BBox3fa& box_o)
    {
      ssef lower(pos_inf);
      ssef upper(neg_inf);
      size_t items = current.size();
      size_t start = current.begin;
      assert(items<=4);
      
      /* allocate leaf node */
      Triangle4i* accel = (Triangle4i*) leafAlloc.malloc(sizeof(Triangle4i));
      *current.parent = bvh->encodeLeaf((char*)accel,listMode ? listMode : 1);

      ssei vgeomID = -1, vprimID = -1;
      Vec3f* v0[4] = { NULL, NULL, NULL, NULL };
      ssei v1 = zero, v2 = zero;
      
      for (size_t i=0; i<items; i++)
      {
	const size_t index = morton[start+i].index;
        const size_t primID = index & encodeMask; 
        const size_t geomID = this->mesh ? this->mesh->id : (index >> encodeShift); 
        const TriangleMesh* mesh = scene->getTriangleMesh(geomID);
	const TriangleMesh::Triangle& tri = mesh->triangle(primID);
	const Vec3fa& p0 = mesh->vertex(tri.v[0]);
        const Vec3fa& p1 = mesh->vertex(tri.v[1]);
        const Vec3fa& p2 = mesh->vertex(tri.v[2]);
        lower = min(lower,(ssef)p0,(ssef)p1,(ssef)p2);
        upper = max(upper,(ssef)p0,(ssef)p1,(ssef)p2);
	vgeomID[i] = geomID;
	vprimID[i] = primID;
	v0[i] = (Vec3f*) mesh->vertexPtr(tri.v[0]); 
	v1[i] = (int*)   mesh->vertexPtr(tri.v[1])-(int*)v0[i]; 
	v2[i] = (int*)   mesh->vertexPtr(tri.v[2])-(int*)v0[i]; 
      }

      for (size_t i=items; i<4; i++)
      {
	vgeomID[i] = -1;
	vprimID[i] = -1;
	v0[i] = v0[0];
	v1[i] = 0; 
	v2[i] = 0;
      }
    
      new (accel) Triangle4i(v0,v1,v2,vgeomID,vprimID,listMode);
      box_o = BBox3fa((Vec3fa)lower,(Vec3fa)upper);
    }
    
    void BVH4BuilderMorton::splitFallback(BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild) const
    {
      const unsigned int center = (current.begin + current.end)/2;
      leftChild.init(current.begin,center);
      rightChild.init(center,current.end);
    }
    
    BBox3fa BVH4BuilderMorton::createLeaf(BuildRecord& current, Allocator& nodeAlloc, Allocator& leafAlloc, size_t threadID)
    {
#if defined(DEBUG)
      if (current.depth > BVH4::maxBuildDepthLeaf) 
        THROW_RUNTIME_ERROR("ERROR: depth limit reached");
#endif
      
      /* create leaf for few primitives */
      if (current.size() <= minLeafSize) {
        BBox3fa bounds;
        createSmallLeaf(current,leafAlloc,threadID,bounds);
        return bounds;
      }
      
      /* first split level */
      BuildRecord record0, record1;
      splitFallback(current,record0,record1);
      
      /* second split level */
      BuildRecord children[4];
      splitFallback(record0,children[0],children[1]);
      splitFallback(record1,children[2],children[3]);
      
      /* allocate node */
      Node* node = (Node*) nodeAlloc.malloc(sizeof(Node)); node->clear();
      *current.parent = bvh->encodeNode(node);
      
      /* recurse into each child */
      BBox3fa bounds0 = empty;
      for (size_t i=0; i<4; i++) {
        children[i].parent = &node->child(i);
        children[i].depth = current.depth+1;
        BBox3fa bounds = createLeaf(children[i],nodeAlloc,leafAlloc,threadID);
        bounds0.extend(bounds);
        node->set(i,bounds);
      }
      BVH4::compact(node); // move empty nodes to the end
      return bounds0;
    }  
    
    __forceinline void BVH4BuilderMorton::split(BuildRecord& current,
                                                BuildRecord& left,
                                                BuildRecord& right) const
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
    
    BBox3fa BVH4BuilderMorton::recurse(BuildRecord& current, Allocator& nodeAlloc, Allocator& leafAlloc, const size_t mode, const size_t threadID) 
    {
      /* stop toplevel recursion at some number of items */
      if (mode == CREATE_TOP_LEVEL && current.size() <= topLevelItemThreshold) {
	state->buildRecords.push_back(current);
        return empty;
      }
      
      __aligned(64) BuildRecord children[BVH4::N];
      
      /* create leaf node */
      if (unlikely(current.depth >= BVH4::maxBuildDepth || current.size() <= minLeafSize)) {
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
          if (children[i].size() <= minLeafSize)
            continue;
          
          /* remember child with largest area */
          if (children[i].size() > bestItems) { 
            bestItems = children[i].size();
            bestChild = i;
          }
        }
        if (bestChild == -1) break;
        
        /*! split best child into left and right child */
        __aligned(64) BuildRecord left, right;
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
        BBox3fa bounds; createSmallLeaf(current,leafAlloc,threadID,bounds); return bounds;
      }
      
      /* allocate node */
      Node* node = (Node*) nodeAlloc.malloc(sizeof(Node)); node->clear();
      *current.parent = bvh->encodeNode(node);
      
      /* recurse into each child */
      BBox3fa bounds0 = empty;
      for (size_t i=0; i<numChildren; i++) 
      {
        children[i].parent = &node->child(i);
        
        if (children[i].size() <= minLeafSize) {
          const BBox3fa bounds = createLeaf(children[i],nodeAlloc,leafAlloc,threadID);
          bounds0.extend(bounds);
          node->set(i,bounds);
        } else {
          const BBox3fa bounds = recurse(children[i],nodeAlloc,leafAlloc,mode,threadID);
          bounds0.extend(bounds);
          node->set(i,bounds);
        }
      }
      return bounds0;
    }
    
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    

    __forceinline BBox3fa BVH4Triangle1BuilderMorton::leafBounds(NodeRef& ref) const
    {
      BBox3fa bounds = empty;
      size_t num; Triangle1* tri = (Triangle1*) ref.leaf(num);
      if (listMode) {
	do {
	  bounds.extend(tri->bounds());
	} while (!((tri++)->last()));
      } else {
	for (size_t i=0; i<num; i++)  // FIXME: have to iterate list also
	  bounds.extend(tri[i].bounds());
      }
      return bounds;
    }
    
    __forceinline BBox3fa BVH4Triangle4BuilderMorton::leafBounds(NodeRef& ref) const
    {
      BBox3fa bounds = empty;
      size_t num; Triangle4* tri = (Triangle4*) ref.leaf(num);
      if (listMode) {
	do {
	  bounds.extend(tri->bounds());
	} while (!((tri++)->last()));
      }
      else
      {
	for (size_t i=0; i<num; i++) 
	  bounds.extend(tri[i].bounds());
      }
      return bounds;
    }

#if defined(__AVX__)
    __forceinline BBox3fa BVH4Triangle8BuilderMorton::leafBounds(NodeRef& ref) const
    {
      BBox3fa bounds = empty;
      size_t num; Triangle8* tri = (Triangle8*) ref.leaf(num);
      if (listMode) {
	do {
	  bounds.extend(tri->bounds());
	} while (!((tri++)->last()));
      }
      else
      {
	for (size_t i=0; i<num; i++) 
	  bounds.extend(tri[i].bounds());
      }
      return bounds;
    }
#endif
    
    __forceinline BBox3fa BVH4Triangle1vBuilderMorton::leafBounds(NodeRef& ref) const
    {
      BBox3fa bounds = empty;
      size_t num; Triangle1v* tri = (Triangle1v*) ref.leaf(num);
      if (listMode) {
	do {
	  bounds.extend(tri->bounds());
	} while (!((tri++)->last()));
      }
      else
      {
	for (size_t i=0; i<num; i++) 
	  bounds.extend(tri[i].bounds());
      }
      return bounds;
    }
    
    __forceinline BBox3fa BVH4Triangle4vBuilderMorton::leafBounds(NodeRef& ref) const
    {
      BBox3fa bounds = empty;
      size_t num; Triangle4v* tri = (Triangle4v*) ref.leaf(num);
      if (listMode) {
	do {
	  bounds.extend(tri->bounds());
	} while (!((tri++)->last()));
      }
      else
      {
	for (size_t i=0; i<num; i++) 
	  bounds.extend(tri[i].bounds());
      }
      return bounds;
    }

    __forceinline BBox3fa BVH4Triangle4iBuilderMorton::leafBounds(NodeRef& ref) const
    {
      BBox3fa bounds = empty;
      size_t num; Triangle4i* tri = (Triangle4i*) ref.leaf(num);
      if (listMode) {
	do {
	  bounds.extend(tri->bounds());
	} while (!((tri++)->last()));
      }
      else
      {
	for (size_t i=0; i<num; i++) 
	  bounds.extend(tri[i].bounds());
      }
      return bounds;
    }
    
    __forceinline BBox3fa BVH4BuilderMorton::nodeBounds(NodeRef& ref) const
    {
      if (ref.isNode())
        return ref.node()->bounds();
      else
        return leafBounds(ref);
    }
    
    BBox3fa BVH4BuilderMorton::refitTopLevel(NodeRef& ref) const
    { 
      /* stop here if we encounter a barrier */
      if (unlikely(ref.isBarrier())) {
        ref.clearBarrier();
        return nodeBounds(ref);
      }
      
      /* return point bound for empty nodes */
      if (unlikely(ref == BVH4::emptyNode))
        return BBox3fa(empty);
      
      /* this is a leaf node */
      if (unlikely(ref.isLeaf()))
	return leafBounds(ref);
      
      /* recurse if this is an internal node */
      Node* node = ref.node();
      const BBox3fa bounds0 = refitTopLevel(node->child(0));
      const BBox3fa bounds1 = refitTopLevel(node->child(1));
      const BBox3fa bounds2 = refitTopLevel(node->child(2));
      const BBox3fa bounds3 = refitTopLevel(node->child(3));
      
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
      return BBox3fa(Vec3fa(lower_x,lower_y,lower_z),
		     Vec3fa(upper_x,upper_y,upper_z));
    }

    void BVH4BuilderMorton::build_sequential_morton(size_t threadIndex, size_t threadCount) 
    {
      /* start measurement */
      double t0 = 0.0f;
      if (g_verbose >= 2) t0 = getSeconds();

      /* compute scene bounds */
      global_bounds = computeBounds();
      bvh->bounds = global_bounds.geomBounds;

      /* compute morton codes */
      size_t dst = 0;
      if (mesh) computeMortonCodes(0,numPrimitives,dst,mesh->id,0,morton);
      else      computeMortonCodes(0,numPrimitives,dst,0,0,morton);
      numPrimitives = dst;

      /* sort morton codes */
      std::sort(&morton[0],&morton[numPrimitives]); // FIXME: use radix sort
      
#if defined(DEBUG)
      for (size_t i=1; i<numPrimitives; i++)
        assert(morton[i-1].code <= morton[i].code);
#endif	    
      
      BuildRecord br;
      br.init(0,numPrimitives);
      br.parent = &bvh->root;
      br.depth = 1;
      
      /* perform first splits in single threaded mode */
      bvh->alloc.clear();
      __aligned(64) Allocator nodeAlloc(&bvh->alloc);
      __aligned(64) Allocator leafAlloc(&bvh->alloc);
      recurse(br,nodeAlloc,leafAlloc,RECURSE,threadIndex);	    
      _mm_sfence(); // make written leaves globally visible
            
      /* stop measurement */
      if (g_verbose >= 2) dt = getSeconds()-t0;
    }
    
    void BVH4BuilderMorton::build_parallel_morton(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      /* start measurement */
      double t0 = 0.0f;
      if (g_verbose >= 2) t0 = getSeconds();

      /* compute scene bounds */
      global_bounds.reset();
      scheduler->dispatchTask( task_computeBounds, this, threadIndex, threadCount );
      bvh->bounds = global_bounds.geomBounds;

      /* calculate initial destination for each thread */
      for (size_t i=0; i<threadCount; i++)
	state->dest[i] = i*numPrimitives/threadCount;

      /* compute morton codes */
      scheduler->dispatchTask( task_computeMortonCodes, this, threadIndex, threadCount );   

      /* calculate new destinations */
      size_t cnt = 0;
      for (size_t i=0; i<threadCount; i++) {
	size_t n = state->dest[i]; state->dest[i] = cnt; cnt += n;
      }
      
      /* if primitive got filtered out, run again */
      if (cnt < numPrimitives) {
	scheduler->dispatchTask( task_computeMortonCodes, this, threadIndex, threadCount );   
	numPrimitives = cnt;
      }
      
      /* padding */
      MortonID32Bit* __restrict__ const dest = (MortonID32Bit*) bvh->alloc.base();
      for (size_t i=numPrimitives; i<( (numPrimitives+7)&(-8) ); i++) {
        dest[i].code  = 0xffffffff; 
        dest[i].index = 0;
      }

      /* sort morton codes */
      barrier.init(threadCount);
      scheduler->dispatchTask( task_radixsort, this, threadIndex, threadCount );

#if defined(DEBUG)
      for (size_t i=1; i<numPrimitives; i++)
        assert(morton[i-1].code <= morton[i].code);
#endif	    

      /* build and extract top-level tree */
      state->buildRecords.clear();
      topLevelItemThreshold = (numPrimitives + threadCount-1)/(2*threadCount);
      
      BuildRecord br;
      br.init(0,numPrimitives);
      br.parent = &bvh->root;
      br.depth = 1;
      
      /* perform first splits in single threaded mode */
      bvh->alloc.clear();
      __aligned(64) Allocator nodeAlloc(&bvh->alloc);
      __aligned(64) Allocator leafAlloc(&bvh->alloc);
      recurse(br,nodeAlloc,leafAlloc,CREATE_TOP_LEVEL,threadIndex);	    
      _mm_sfence(); // make written leaves globally visible

      /* sort all subtasks by size */
      std::sort(state->buildRecords.begin(),state->buildRecords.end(),BuildRecord::Greater());

      /* build sub-trees */
      state->taskCounter = 0;
      //state->workStack.reset();
      scheduler->dispatchTask( task_recurseSubMortonTrees, this, threadIndex, threadCount );
      
      /* refit toplevel part of tree */
      refitTopLevel(bvh->root);

      /* stop measurement */
      if (g_verbose >= 2) dt = getSeconds()-t0;
    }

    Builder* BVH4Triangle1BuilderMorton  (void* bvh, Scene* scene, size_t mode) { return new class BVH4Triangle1BuilderMorton ((BVH4*)bvh,scene,mode); }
    Builder* BVH4Triangle4BuilderMorton  (void* bvh, Scene* scene, size_t mode) { return new class BVH4Triangle4BuilderMorton ((BVH4*)bvh,scene,mode); }
#if defined(__AVX__)
    Builder* BVH4Triangle8BuilderMorton  (void* bvh, Scene* scene, size_t mode) { return new class BVH4Triangle8BuilderMorton ((BVH4*)bvh,scene,mode); }
#endif
    Builder* BVH4Triangle1vBuilderMorton (void* bvh, Scene* scene, size_t mode) { return new class BVH4Triangle1vBuilderMorton((BVH4*)bvh,scene,mode); }
    Builder* BVH4Triangle4vBuilderMorton (void* bvh, Scene* scene, size_t mode) { return new class BVH4Triangle4vBuilderMorton((BVH4*)bvh,scene,mode); }
    Builder* BVH4Triangle4iBuilderMorton (void* bvh, Scene* scene, size_t mode) { return new class BVH4Triangle4iBuilderMorton((BVH4*)bvh,scene,mode); }

    Builder* BVH4Triangle1MeshBuilderMorton  (void* bvh, TriangleMesh* mesh, size_t mode) { return new class BVH4Triangle1BuilderMorton ((BVH4*)bvh,mesh,mode); }
    Builder* BVH4Triangle4MeshBuilderMorton  (void* bvh, TriangleMesh* mesh, size_t mode) { return new class BVH4Triangle4BuilderMorton ((BVH4*)bvh,mesh,mode); }
#if defined(__AVX__)
    Builder* BVH4Triangle8MeshBuilderMorton  (void* bvh, TriangleMesh* mesh, size_t mode) { return new class BVH4Triangle8BuilderMorton ((BVH4*)bvh,mesh,mode); }
#endif
    Builder* BVH4Triangle1vMeshBuilderMorton (void* bvh, TriangleMesh* mesh, size_t mode) { return new class BVH4Triangle1vBuilderMorton((BVH4*)bvh,mesh,mode); }
    Builder* BVH4Triangle4vMeshBuilderMorton (void* bvh, TriangleMesh* mesh, size_t mode) { return new class BVH4Triangle4vBuilderMorton((BVH4*)bvh,mesh,mode); }
    Builder* BVH4Triangle4iMeshBuilderMorton (void* bvh, TriangleMesh* mesh, size_t mode) { return new class BVH4Triangle4iBuilderMorton((BVH4*)bvh,mesh,mode); }
  }
}

  
  
