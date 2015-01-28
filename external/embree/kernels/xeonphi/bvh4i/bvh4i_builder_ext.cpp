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

#include "bvh4i/bvh4i.h"
#include "bvh4i/bvh4i_builder.h"

#define PRESPLIT_SPACE_FACTOR         0.30f
#define PRESPLIT_AREA_THRESHOLD      20.0f
#define PRESPLIT_MIN_AREA             0.01f
#define NUM_PRESPLITS_PER_TRIANGLE    16
#define PRESPLITS_TREE_DEPTH          4
#define NUM_PRESPLIT_IDS_PER_BLOCK    8

#define DBG(x) 
#define TIMER(x) 

#define L1_PREFETCH_ITEMS 2
#define L2_PREFETCH_ITEMS 16

namespace embree
{

  /* =================================================================================== */
  /* =================================================================================== */
  /* =================================================================================== */

  struct __aligned(16) PreSplitID
  {
    unsigned int code;
    float sah;
    unsigned int groupID;
    unsigned int primID;
            
    __forceinline unsigned int get(const unsigned int shift, const unsigned and_mask) const {
      return (code >> shift) & and_mask;
    }

    __forceinline unsigned int getByte(const size_t b) const {
      assert(b < 4);
      const unsigned char *__restrict const ptr = (const unsigned char*)&code;
      return ptr[b];
    }
            
    __forceinline friend std::ostream &operator<<(std::ostream &o, const PreSplitID& pid) {
      o << "code " << pid.code << " boxSAH = " << pid.sah << " groupID " << pid.groupID << " primID " << pid.primID;
      return o;
    }
    __forceinline bool operator<(const PreSplitID &pid) const { return code < pid.code; } 
    __forceinline bool operator>(const PreSplitID &pid) const { return code > pid.code; } 
  };

  void BVH4iBuilderPreSplits::printBuilderName()
  {
    std::cout << "building BVH4i with presplits-based SAH builder (MIC) ... " << std::endl;    
  }

  void BVH4iBuilderPreSplits::allocateData(const size_t threadCount, const size_t totalNumPrimitives)
  {
    DBG(PING);
    size_t numPrimitivesOld = numPrimitives;
    numPrimitives = totalNumPrimitives;
    DBG(DBG_PRINT(numPrimitives));

    if (numPrimitivesOld != numPrimitives)
      {
	const size_t preSplitPrims = (size_t)((float)numPrimitives * PRESPLIT_SPACE_FACTOR);
	const size_t numPrims = numPrimitives+preSplitPrims;
	const size_t minAllocNodes = numPrims ? threadCount * ALLOCATOR_NODE_BLOCK_SIZE * 4: 16;
	const size_t numNodes = max((size_t)(numPrims * BVH_NODE_PREALLOC_FACTOR),minAllocNodes);

	numMaxPrimitives = numPrims;
	numMaxPreSplits  = numPrims - numPrimitives;

	if (g_verbose >= 2)
	  {
	    DBG_PRINT(numPrimitives);
	    DBG_PRINT(numMaxPrimitives);
	    DBG_PRINT(numMaxPreSplits);
	  };

	allocateMemoryPools(numPrims,numNodes);
      }    
  }


  __forceinline void splitTri(const PrimRef& prim, int dim, float pos, const Vec3fa& a, const Vec3fa& b, const Vec3fa& c, PrimRef& left_o, PrimRef& right_o)
  {
    BBox3f left = empty, right = empty;
    const Vec3fa v[3] = { a,b,c };

    /* clip triangle to left and right box by processing all edges */
    Vec3fa v1 = v[2];
    for (size_t i=0; i<3; i++)
      {
	Vec3fa v0 = v1; v1 = v[i];
	float v0d = v0[dim], v1d = v1[dim];
      
	if (v0d <= pos) left. extend(v0); // this point is on left side
	if (v0d >= pos) right.extend(v0); // this point is on right side

	if ((v0d < pos && pos < v1d) || (v1d < pos && pos < v0d)) // the edge crosses the splitting location
	  {
	    assert((v1d-v0d) != 0.0f);
	    Vec3fa c = v0 + (pos-v0d)/(v1d-v0d)*(v1-v0);
	    left.extend(c);
	    right.extend(c);
	  }
      }
    assert(!left.empty());  // happens if split does not hit triangle
    assert(!right.empty()); // happens if split does not hit triangle

    /* safe clip against current bounds */
    BBox3f bounds = prim.bounds();
    BBox3f cleft(min(max(left.lower,bounds.lower),bounds.upper),
                 max(min(left.upper,bounds.upper),bounds.lower));
    BBox3f cright(min(max(right.lower,bounds.lower),bounds.upper),
                  max(min(right.upper,bounds.upper),bounds.lower));

    left_o  = PrimRef(cleft, prim.geomID(), prim.primID());
    right_o = PrimRef(cright,prim.geomID(), prim.primID());
  }

  __forceinline mic_f box_sah( const mic_f &b_min,
			       const mic_f &b_max) 
  { 
    const mic_f d = b_max - b_min;
    const mic_f d_x = swAAAA(d);
    const mic_f d_y = swBBBB(d);
    const mic_f d_z = swCCCC(d);
    return (d_x*(d_y+d_z)+d_y*d_z)*2.0f; 
  }

  __forceinline mic_f box_sah( const PrimRef &r) 
  { 
    const mic_f bmin = broadcast4to16f(&r.lower);
    const mic_f bmax = broadcast4to16f(&r.upper);
    return box_sah(bmin,bmax);
  }

  __forceinline mic_f tri_sah( const mic_f &v0,
			       const mic_f &v1,
			       const mic_f &v2) 
  {
    const mic_f n = lcross_xyz(v1-v0,v2-v0);
    return sqrt(ldot3_xyz(n,n)) * 0.5f;
  }
  
  __forceinline size_t getMaxDim(const PrimRef &p)
  {
    Vec3fa diag = p.upper - p.lower;
    size_t index = 0;
    for (size_t i=1;i<3;i++)
      if (diag[i] > diag[index])
	index = i;
    return index;
  }

  void subdivideTriangle(const PrimRef &primBounds,
			 const Vec3fa& vtxA,
			 const Vec3fa& vtxB,
			 const Vec3fa& vtxC,
			 const size_t depth,
			 AlignedAtomicCounter32 &counter,
			 PrimRef *__restrict__ prims)
  {
    DBG(
	DBG_PRINT(primBounds);
	DBG_PRINT(depth);
	);

    if (depth == 0) 
      {	    
	const unsigned int index = counter.inc();
	prims[index] = primBounds;
      }
    else
      {
	const size_t dim = getMaxDim(primBounds);
	
	PrimRef left,right;

	const float pos = (primBounds.upper[dim] + primBounds.lower[dim]) * 0.5f;
	splitTri(primBounds,dim,pos,vtxA,vtxB,vtxC,left,right);

	DBG(
	    DBG_PRINT(left);
	    DBG_PRINT(right);
	    );

	subdivideTriangle( left ,vtxA,vtxB,vtxC, depth-1,counter,prims);
	subdivideTriangle( right,vtxA,vtxB,vtxC, depth-1,counter,prims);
      }
  }


  void BVH4iBuilderPreSplits::computePrimRefs(const size_t threadIndex, const size_t threadCount)
  {
    TIMER(double msec = 0.0);
    TIMER(msec = getSeconds());

    dest0.reset(0);
    dest1.reset(0);

    LockStepTaskScheduler::dispatchTask( task_countAndComputePrimRefsPreSplits, this, threadIndex, threadCount );

    /* === padding to 8-wide blocks === */
    const unsigned int preSplits        = dest1;
    const unsigned int preSplits_padded = ((preSplits+NUM_PRESPLIT_IDS_PER_BLOCK-1)&(-NUM_PRESPLIT_IDS_PER_BLOCK));
    PreSplitID* __restrict__ const dest = (PreSplitID*)accel;
    
    for (size_t i=preSplits; i<preSplits_padded; i++) {
      dest[i].code  = 0xffffffff; 
      dest[i].sah = 0;
      dest[i].groupID = 0;
      dest[i].primID = 0;
    }

    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "task_countAndComputePrimRefsPreSplits " << 1000. * msec << " ms" << std::endl << std::flush);


    const size_t step = (numMaxPreSplits+NUM_PRESPLITS_PER_TRIANGLE-1) / NUM_PRESPLITS_PER_TRIANGLE;
    const size_t startPreSplits = (step <= preSplits) ? preSplits - step : preSplits;


    DBG(
	DBG_PRINT(numPrimitives);
	DBG_PRINT(preSplits);
	DBG_PRINT(preSplits_padded);
	DBG_PRINT(dest0);
	DBG_PRINT(dest1);
	DBG_PRINT(step);
	DBG_PRINT(startPreSplits);
	);

    TIMER(msec = getSeconds());

    LockStepTaskScheduler::dispatchTask( task_radixSortPreSplitIDs, this, threadIndex, threadCount );

    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "task_radixSortPreSplitIDs " << 1000. * msec << " ms" << std::endl << std::flush);


    TIMER(msec = getSeconds());

    LockStepTaskScheduler::dispatchTask( task_computePrimRefsFromPreSplitIDs, this, threadIndex, threadCount );

    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "task_computePrimRefsFromPreSplitIDs " << 1000. * msec << " ms" << std::endl << std::flush);


    numPrimitives = dest0;

    DBG(
	DBG_PRINT(numPrimitives);
	);

  }



  void BVH4iBuilderPreSplits::countAndComputePrimRefsPreSplits(const size_t threadID, const size_t numThreads) 
  {
    DBG(PING);
    const size_t numGroups = scene->size();
    const size_t startID = (threadID+0)*numPrimitives/numThreads;
    const size_t endID   = (threadID+1)*numPrimitives/numThreads;
    

    // === find first group containing startID ===
    unsigned int startGroup=0, numSkipped = 0;
    for (; startGroup<numGroups; startGroup++) {       
      if (unlikely(scene->get(startGroup) == NULL)) continue;
      if (unlikely(scene->get(startGroup)->type != TRIANGLE_MESH)) continue;
      const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(startGroup);
      if (unlikely(!mesh->isEnabled())) continue;
      if (unlikely(mesh->numTimeSteps != 1)) continue;

      const size_t numTriangles = mesh->numTriangles;
      if (numSkipped + numTriangles > startID) break;
      numSkipped += numTriangles;
    }

    // === start with first group containing startID ===
    mic_f bounds_scene_min((float)pos_inf);
    mic_f bounds_scene_max((float)neg_inf);
    mic_f bounds_centroid_min((float)pos_inf);
    mic_f bounds_centroid_max((float)neg_inf);

    unsigned int numTrisNoPreSplit = 0;
    unsigned int numTrisPreSplit = 0;

    // === determine presplit candidates ===
    {
      unsigned int currentID = startID;
      unsigned int offset = startID - numSkipped;
      for (unsigned int g=startGroup; g<numGroups; g++) 
	{
	  if (unlikely(scene->get(g) == NULL)) continue;
	  if (unlikely(scene->get(g)->type != TRIANGLE_MESH)) continue;
	  const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(g);
	  if (unlikely(!mesh->isEnabled())) continue;
	  if (unlikely(mesh->numTimeSteps != 1)) continue;

	  for (unsigned int i=offset; i<mesh->numTriangles && currentID < endID; i++, currentID++)	 
	    { 			    
	      const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(i);
	      prefetch<PFHINT_L2>(&tri + L2_PREFETCH_ITEMS);
	      prefetch<PFHINT_L1>(&tri + L1_PREFETCH_ITEMS);

	      const float *__restrict__ const vptr0 = (float*)&mesh->vertex(tri.v[0]);
	      const float *__restrict__ const vptr1 = (float*)&mesh->vertex(tri.v[1]);
	      const float *__restrict__ const vptr2 = (float*)&mesh->vertex(tri.v[2]);

	      const mic_f v0 = broadcast4to16f(vptr0);
	      const mic_f v1 = broadcast4to16f(vptr1);
	      const mic_f v2 = broadcast4to16f(vptr2);

	      mic_f bmin = min(min(v0,v1),v2);
	      mic_f bmax = max(max(v0,v1),v2);


	      const mic_f area_tri = tri_sah(v0,v1,v2);
	      const mic_f area_box = box_sah(bmin,bmax);
	      const mic_f factor = area_box * rcp(area_tri);

	      DBG(
		  DBG_PRINT(area_tri);
		  DBG_PRINT(area_box);
		  DBG_PRINT(factor);
		  );

	      const mic_m m_factor = factor > PRESPLIT_AREA_THRESHOLD;
	      const mic_m m_sah_zero = area_box > PRESPLIT_MIN_AREA;

	      if (any(m_factor & m_sah_zero)) 
		numTrisPreSplit++;
	      else
		numTrisNoPreSplit++;   

	      bounds_scene_min = min(bounds_scene_min,bmin);
	      bounds_scene_max = max(bounds_scene_max,bmax);
	      const mic_f centroid2 = bmin+bmax;
	      bounds_centroid_min = min(bounds_centroid_min,centroid2);
	      bounds_centroid_max = max(bounds_centroid_max,centroid2);
	    }
	  if (currentID == endID) break;
	  offset = 0;
	}
    }

    /* update global bounds */
    Centroid_Scene_AABB bounds;
    
    store4f(&bounds.centroid2.lower,bounds_centroid_min);
    store4f(&bounds.centroid2.upper,bounds_centroid_max);
    store4f(&bounds.geometry.lower,bounds_scene_min);
    store4f(&bounds.geometry.upper,bounds_scene_max);

    global_bounds.extend_atomic(bounds);    


    const unsigned int no_presplit_index = dest0.add(numTrisNoPreSplit);
    const unsigned int presplit_index    = dest1.add(numTrisPreSplit);

    PrimRef *presplits = (PrimRef*)accel;

    PrimRef    *__restrict__ no_presplit_prims = this->prims + no_presplit_index;
    PreSplitID *__restrict__ presplitIDs       = (PreSplitID*)accel + presplit_index;


    // === put triangles into prepsplit/no-presplit bucket ===
    {
      unsigned int currentID = startID;
      unsigned int offset = startID - numSkipped;
      for (unsigned int g=startGroup; g<numGroups; g++) 
	{
	  if (unlikely(scene->get(g) == NULL)) continue;
	  if (unlikely(scene->get(g)->type != TRIANGLE_MESH)) continue;
	  const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(g);
	  if (unlikely(!mesh->isEnabled())) continue;
	  if (unlikely(mesh->numTimeSteps != 1)) continue;

	  for (unsigned int i=offset; i<mesh->numTriangles && currentID < endID; i++, currentID++)	 
	    { 			    
	      const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(i);
	      prefetch<PFHINT_L2>(&tri + L2_PREFETCH_ITEMS);
	      prefetch<PFHINT_L1>(&tri + L1_PREFETCH_ITEMS);

	      const float *__restrict__ const vptr0 = (float*)&mesh->vertex(tri.v[0]);
	      const float *__restrict__ const vptr1 = (float*)&mesh->vertex(tri.v[1]);
	      const float *__restrict__ const vptr2 = (float*)&mesh->vertex(tri.v[2]);

	      const mic_f v0 = broadcast4to16f(vptr0);
	      const mic_f v1 = broadcast4to16f(vptr1);
	      const mic_f v2 = broadcast4to16f(vptr2);

	      mic_f bmin = min(min(v0,v1),v2);
	      mic_f bmax = max(max(v0,v1),v2);


	      const mic_f area_tri = tri_sah(v0,v1,v2);
	      const mic_f area_box = box_sah(bmin,bmax);
	      const mic_f factor = area_box * rcp(area_tri);

	      DBG(
		  DBG_PRINT(area_tri);
		  DBG_PRINT(area_box);
		  DBG_PRINT(factor);
		  );

	      const mic_m m_factor = factor > PRESPLIT_AREA_THRESHOLD;
	      const mic_m m_sah_zero = area_box > PRESPLIT_MIN_AREA;

	      if (any(m_factor & m_sah_zero)) 
		{
		  prefetch<PFHINT_L2EX>(presplitIDs + 4*4);

		  store1f(&presplitIDs->code,area_box);
		  store1f(&presplitIDs->sah,factor);
		  presplitIDs->groupID = g;
		  presplitIDs->primID  = i;
		  presplitIDs++;
		}
	      else
		{
		  prefetch<PFHINT_L2EX>(no_presplit_prims + L2_PREFETCH_ITEMS);

		  store4f(&no_presplit_prims->lower,bmin);
		  store4f(&no_presplit_prims->upper,bmax);	
		  no_presplit_prims->lower.a = g; 
		  no_presplit_prims->upper.a = i;
		  no_presplit_prims++;

		}

	    }
	  if (currentID == endID) break;
	  offset = 0;
	}
    }
  }


  void BVH4iBuilderPreSplits::radixSortPreSplitIDs(const size_t threadID, const size_t numThreads)
  {
    const unsigned int preSplits = dest1;
    const size_t numBlocks = (preSplits+NUM_PRESPLIT_IDS_PER_BLOCK-1) / NUM_PRESPLIT_IDS_PER_BLOCK;
    const size_t startID   = ((threadID+0)*numBlocks/numThreads) * NUM_PRESPLIT_IDS_PER_BLOCK;
    const size_t endID     = ((threadID+1)*numBlocks/numThreads) * NUM_PRESPLIT_IDS_PER_BLOCK;
    assert(startID % NUM_PRESPLIT_IDS_PER_BLOCK == 0);
    assert(endID % NUM_PRESPLIT_IDS_PER_BLOCK == 0);

    assert(((numThreads)*numBlocks/numThreads) * NUM_PRESPLIT_IDS_PER_BLOCK == ((preSplits+7)&(-8)));


    PreSplitID* __restrict__ presplitID[2];
    presplitID[0] = (PreSplitID*)accel; 
    presplitID[1] = (PreSplitID*)accel + ((preSplits+NUM_PRESPLIT_IDS_PER_BLOCK-1)&(-NUM_PRESPLIT_IDS_PER_BLOCK));


    /* we need 4 iterations to process all 32 bits */
    for (size_t b=0; b<4; b++)
      {
	const PreSplitID* __restrict__ const src = (PreSplitID*)presplitID[((b+0)%2)];
	PreSplitID*       __restrict__ const dst = (PreSplitID*)presplitID[((b+1)%2)];

	__assume_aligned(&radixCount[threadID][0],64);
      
	/* count how many items go into the buckets */

#pragma unroll(16)
	for (size_t i=0; i<16; i++)
	  store16i(&radixCount[threadID][i*16],mic_i::zero());


	for (size_t i=startID; i<endID; i+=NUM_PRESPLIT_IDS_PER_BLOCK) {
	  prefetch<PFHINT_NT>(&src[i+L1_PREFETCH_ITEMS]);
	  prefetch<PFHINT_L2>(&src[i+L2_PREFETCH_ITEMS]);
	
#pragma unroll(NUM_PRESPLIT_IDS_PER_BLOCK)
	  for (unsigned long j=0;j<NUM_PRESPLIT_IDS_PER_BLOCK;j++)
	    {
	      const unsigned int index = src[i+j].getByte(b);
	      radixCount[threadID][index]++;
	    }
	}

	LockStepTaskScheduler::syncThreads(threadID,numThreads);


	/* calculate total number of items for each bucket */


	mic_i count[16];
#pragma unroll(16)
	for (size_t i=0; i<16; i++)
	  count[i] = mic_i::zero();


	for (size_t i=0; i<threadID; i++)
#pragma unroll(16)
	  for (size_t j=0; j<16; j++)
	    count[j] += load16i((int*)&radixCount[i][j*16]);
      
	__aligned(64) unsigned int inner_offset[RADIX_BUCKETS];

#pragma unroll(16)
	for (size_t i=0; i<16; i++)
	  store16i(&inner_offset[i*16],count[i]);

#pragma unroll(16)
	for (size_t i=0; i<16; i++)
	  count[i] = load16i((int*)&inner_offset[i*16]);

	for (size_t i=threadID; i<numThreads; i++)
#pragma unroll(16)
	  for (size_t j=0; j<16; j++)
	    count[j] += load16i((int*)&radixCount[i][j*16]);	  

	__aligned(64) unsigned int total[RADIX_BUCKETS];

#pragma unroll(16)
	for (size_t i=0; i<16; i++)
	  store16i(&total[i*16],count[i]);

	__aligned(64) unsigned int offset[RADIX_BUCKETS];

	/* calculate start offset of each bucket */
	offset[0] = 0;
	for (size_t i=1; i<RADIX_BUCKETS; i++)    
	  offset[i] = offset[i-1] + total[i-1];
      
	/* calculate start offset of each bucket for this thread */

#pragma unroll(RADIX_BUCKETS)
	for (size_t j=0; j<RADIX_BUCKETS; j++)
          offset[j] += inner_offset[j];

	/* copy items into their buckets */
	for (size_t i=startID; i<endID; i+=NUM_PRESPLIT_IDS_PER_BLOCK) {
	  prefetch<PFHINT_NT>(&src[i+L1_PREFETCH_ITEMS]);
	  prefetch<PFHINT_L2>(&src[i+L2_PREFETCH_ITEMS]);

#pragma nounroll
	  for (unsigned long j=0;j<NUM_PRESPLIT_IDS_PER_BLOCK;j++)
	    {
	      const unsigned int index = src[i+j].getByte(b);
	      assert(index < RADIX_BUCKETS);
	      dst[offset[index]] = src[i+j];
	      prefetch<PFHINT_L2EX>(&dst[offset[index]+L1_PREFETCH_ITEMS]);
	      offset[index]++;
	    }
	  evictL2(&src[i]);
	}

	if (b<3) LockStepTaskScheduler::syncThreads(threadID,numThreads);

      }
  }

  void BVH4iBuilderPreSplits::computePrimRefsFromPreSplitIDs(const size_t threadID, const size_t numThreads) 
  {
    const size_t preSplitIDs    = dest1;
    const size_t step           = (numMaxPreSplits+NUM_PRESPLITS_PER_TRIANGLE-1) / NUM_PRESPLITS_PER_TRIANGLE;
    const size_t startPreSplits = (step <= preSplitIDs) ? preSplitIDs - step : preSplitIDs;

    // === no pre-splits ==
    {
      const size_t startID = (threadID+0)*startPreSplits/numThreads;
      const size_t endID   = (threadID+1)*startPreSplits/numThreads;
      const size_t items   = endID - startID;
      const unsigned int d_index = dest0.add(items);

      PrimRef    *__restrict__ prims = this->prims + d_index;
      PreSplitID *__restrict__ presplitIDs = (PreSplitID*)accel;
      
      for (size_t i=startID;i<endID;i++)
	{
	  prefetch<PFHINT_L2>(&presplitIDs[i] + L2_PREFETCH_ITEMS);
	  prefetch<PFHINT_L1>(&presplitIDs[i] + L1_PREFETCH_ITEMS);

	  const size_t geomID = presplitIDs[i].groupID;
	  const size_t primID = presplitIDs[i].primID;

	  const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
	  const TriangleMeshScene::TriangleMesh::Triangle & tri = mesh->triangle(primID);

	  const float *__restrict__ const vptr0 = (float*)&mesh->vertex(tri.v[0]);
	  const float *__restrict__ const vptr1 = (float*)&mesh->vertex(tri.v[1]);
	  const float *__restrict__ const vptr2 = (float*)&mesh->vertex(tri.v[2]);

	  const mic_f v0 = broadcast4to16f(vptr0);
	  const mic_f v1 = broadcast4to16f(vptr1);
	  const mic_f v2 = broadcast4to16f(vptr2);

	  prefetch<PFHINT_L1EX>(prims);

	  mic_f bmin = min(min(v0,v1),v2);
	  mic_f bmax = max(max(v0,v1),v2);

	  prefetch<PFHINT_L2EX>(prims + L2_PREFETCH_ITEMS);

	  store4f(&prims->lower,bmin);
	  store4f(&prims->upper,bmax);	
	  prims->lower.a = geomID; 
	  prims->upper.a = primID;
	  prims++;
	}
    }

    // === performing pre-splitting ==
    {
      AlignedAtomicCounter32 counter;
      counter.reset(0);

      const unsigned int numPreSplits = preSplitIDs - startPreSplits;
      const size_t startID = (threadID+0)*numPreSplits/numThreads;
      const size_t endID   = (threadID+1)*numPreSplits/numThreads;
      const size_t items   = (endID - startID)*NUM_PRESPLITS_PER_TRIANGLE;
      const unsigned int d_index = dest0.add(items);

      PrimRef    *__restrict__ prims = this->prims + d_index;
      PreSplitID *__restrict__ presplitIDs = (PreSplitID*)accel + startPreSplits;
      
      for (size_t i=startID;i<endID;i++)
	{
	  prefetch<PFHINT_L2>(&presplitIDs[i] + L2_PREFETCH_ITEMS);
	  prefetch<PFHINT_L1>(&presplitIDs[i] + L1_PREFETCH_ITEMS);

	  const size_t geomID = presplitIDs[i].groupID;
	  const size_t primID = presplitIDs[i].primID;

	  const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
	  const TriangleMeshScene::TriangleMesh::Triangle & tri = mesh->triangle(primID);

	  const float *__restrict__ const vptr0 = (float*)&mesh->vertex(tri.v[0]);
	  const float *__restrict__ const vptr1 = (float*)&mesh->vertex(tri.v[1]);
	  const float *__restrict__ const vptr2 = (float*)&mesh->vertex(tri.v[2]);

	  Vec3fa vtxA = *(Vec3fa*)vptr0;
	  Vec3fa vtxB = *(Vec3fa*)vptr1;
	  Vec3fa vtxC = *(Vec3fa*)vptr2;

	  BBox3f bounds = empty;
	  bounds.extend(vtxA);
	  bounds.extend(vtxB);
	  bounds.extend(vtxC);
	  PrimRef primBounds = PrimRef(bounds,geomID,primID);

	  subdivideTriangle(primBounds,vtxA,vtxB,vtxC,PRESPLITS_TREE_DEPTH,counter,prims);
	}
      assert(counter == items);
    }

  }


  /* =================================================================================== */
  /* =================================================================================== */
  /* =================================================================================== */

  void BVH4iBuilderVirtualGeometry::printBuilderName()
  {
    std::cout << "building BVH4i with Virtual Geometry SAH builder (MIC) ... " << std::endl;    
  }

  size_t BVH4iBuilderVirtualGeometry::getNumPrimitives()
  {
    /* count total number of virtual objects */
    size_t numVirtualObjects = 0;       
    for (size_t i=0;i<scene->size();i++)
      {
	if (unlikely(scene->get(i) == NULL)) continue;
	if (unlikely((scene->get(i)->type != USER_GEOMETRY) && (scene->get(i)->type != INSTANCES))) continue;
	if (unlikely(!scene->get(i)->isEnabled())) continue;
        UserGeometryScene::Base* geom = (UserGeometryScene::Base*) scene->get(i);
	numVirtualObjects += geom->size();
      }
    return numVirtualObjects;	
  }

  void BVH4iBuilderVirtualGeometry::computePrimRefs(const size_t threadIndex, const size_t threadCount)
  {
    DBG(PING);
    LockStepTaskScheduler::dispatchTask( task_computePrimRefsVirtualGeometry, this, threadIndex, threadCount );	
  }

  void BVH4iBuilderVirtualGeometry::createAccel(const size_t threadIndex, const size_t threadCount)
  {
    DBG(PING);
    LockStepTaskScheduler::dispatchTask( task_createVirtualGeometryAccel, this, threadIndex, threadCount );
  }

  void BVH4iBuilderVirtualGeometry::computePrimRefsVirtualGeometry(const size_t threadID, const size_t numThreads) 
  {
    DBG(PING);

    const size_t numTotalGroups = scene->size();

    /* count total number of virtual objects */
    const size_t numVirtualObjects = numPrimitives;
    const size_t startID   = (threadID+0)*numVirtualObjects/numThreads;
    const size_t endID     = (threadID+1)*numVirtualObjects/numThreads; 

    DBG(
	DBG_PRINT(numTotalGroups);
	DBG_PRINT(numVirtualObjects);
	DBG_PRINT(startID);
	DBG_PRINT(endID);
	);
    
    PrimRef *__restrict__ const prims     = this->prims;

    // === find first group containing startID ===
    unsigned int g=0, numSkipped = 0;
    for (; g<numTotalGroups; g++) {       
      if (unlikely(scene->get(g) == NULL)) continue;
      if (unlikely((scene->get(g)->type != USER_GEOMETRY) && (scene->get(g)->type != INSTANCES))) continue;
      if (unlikely(!scene->get(g)->isEnabled())) continue;
      const UserGeometryScene::Base* const geom = (UserGeometryScene::Base*) scene->get(g);
      const size_t numPrims = geom->size();
      if (numSkipped + numPrims > startID) break;
      numSkipped += numPrims;
    }

    /* start with first group containing startID */
    mic_f bounds_scene_min((float)pos_inf);
    mic_f bounds_scene_max((float)neg_inf);
    mic_f bounds_centroid_min((float)pos_inf);
    mic_f bounds_centroid_max((float)neg_inf);

    unsigned int num = 0;
    unsigned int currentID = startID;
    unsigned int offset = startID - numSkipped;

    for (; g<numTotalGroups; g++) 
      {
	if (unlikely(scene->get(g) == NULL)) continue;
	if (unlikely((scene->get(g)->type != USER_GEOMETRY ) && (scene->get(g)->type != INSTANCES))) continue;
	if (unlikely(!scene->get(g)->isEnabled())) continue;

	UserGeometryScene::Base *virtual_geometry = (UserGeometryScene::Base *)scene->get(g);

        size_t N = virtual_geometry->size();
        for (unsigned int i=offset; i<N && currentID < endID; i++, currentID++)	 
	  { 			    
	    const BBox3f bounds = virtual_geometry->bounds(i);
	    const mic_f bmin = broadcast4to16f(&bounds.lower); 
	    const mic_f bmax = broadcast4to16f(&bounds.upper);

	    DBG(
		DBG_PRINT(currentID);
		DBG_PRINT(bmin);
		DBG_PRINT(bmax);
		);
          
	    bounds_scene_min = min(bounds_scene_min,bmin);
	    bounds_scene_max = max(bounds_scene_max,bmax);
	    const mic_f centroid2 = bmin+bmax;
	    bounds_centroid_min = min(bounds_centroid_min,centroid2);
	    bounds_centroid_max = max(bounds_centroid_max,centroid2);

	    store4f(&prims[currentID].lower,bmin);
	    store4f(&prims[currentID].upper,bmax);	
	    prims[currentID].lower.a = g;
	    prims[currentID].upper.a = i;
	  }
        if (currentID == endID) break;
        offset = 0;
      }

    /* update global bounds */
    Centroid_Scene_AABB bounds;
    
    store4f(&bounds.centroid2.lower,bounds_centroid_min);
    store4f(&bounds.centroid2.upper,bounds_centroid_max);
    store4f(&bounds.geometry.lower,bounds_scene_min);
    store4f(&bounds.geometry.upper,bounds_scene_max);

    global_bounds.extend_atomic(bounds);    
  }


  void BVH4iBuilderVirtualGeometry::createVirtualGeometryAccel(const size_t threadID, const size_t numThreads)
  {
    DBG(PING);

    const size_t startID = (threadID+0)*numPrimitives/numThreads;
    const size_t endID   = (threadID+1)*numPrimitives/numThreads;

    AccelSetItem *acc = (AccelSetItem*)accel + startID;

    const PrimRef* __restrict__  bptr = prims + startID;

    for (size_t j=startID; j<endID; j++, bptr++, acc++)
      {
	prefetch<PFHINT_NT>(bptr + L1_PREFETCH_ITEMS);
	prefetch<PFHINT_L2>(bptr + L2_PREFETCH_ITEMS);
	assert(bptr->geomID() < scene->size() );
        AccelSet* _accel = (AccelSet*)(UserGeometryScene::Base *) scene->get( bptr->geomID() );
	acc->accel = _accel;
        acc->item = bptr->primID();
      }
  }

};
