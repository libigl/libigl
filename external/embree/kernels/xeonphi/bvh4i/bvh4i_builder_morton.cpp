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

#include "bvh4i_builder_morton.h"
#include "builders/builder_util.h"
#include "bvh4i_rotate.h"

#define MORTON_BVH4I_NODE_PREALLOC_FACTOR   0.9f
#define NUM_MORTON_IDS_PER_BLOCK            8
#define SINGLE_THREADED_BUILD_THRESHOLD     (MAX_MIC_THREADS*64)

//#define PROFILE
#define PROFILE_ITERATIONS 200

#define TIMER(x) 
#define DBG(x) 

#define L1_PREFETCH_ITEMS 8
#define L2_PREFETCH_ITEMS 44

namespace embree 
{
#if defined(DEBUG)
  extern AtomicMutex mtx;
#endif


  // =======================================================================================================
  // =======================================================================================================
  // =======================================================================================================

  __aligned(64) static double dt = 0.0f;

  BVH4iBuilderMorton::BVH4iBuilderMorton (BVH4i* bvh, void* geometry, const bool tree_rotations)
  : bvh(bvh), scene((Scene*)geometry), topLevelItemThreshold(0), encodeShift(0), encodeMask(0), numBuildRecords(0), 
    morton(NULL), node(NULL), accel(NULL), numGroups(0), numPrimitives(0), numNodes(0), numAllocatedNodes(0), size_morton(0), size_node(0), size_accel(0), numPrimitivesOld(-1), enableTreeRotations(tree_rotations)
  {
  }

  BVH4iBuilderMorton::~BVH4iBuilderMorton()
  {
    if (morton) {
      assert(size_morton > 0);
      os_free(morton,size_morton);
    }
  }


  void BVH4iBuilderMorton::initEncodingAllocateData()
  {
    /* calculate total number of primrefs */
    numGroups     = scene->size();
    numPrimitives = 0;


    size_t maxPrimsPerGroup = 0;
    for (size_t group=0; group<numGroups; group++) 
      {
	if (unlikely(scene->get(group) == NULL)) continue;
	if (scene->get(group)->type != TRIANGLE_MESH) continue;
	const TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(group);
	if (unlikely(!mesh->isEnabled())) continue;
	if (unlikely(mesh->numTimeSteps != 1)) continue;

	maxPrimsPerGroup = max(maxPrimsPerGroup,mesh->numTriangles);
	numPrimitives   += mesh->numTriangles;
      }

    /* calculate groupID, primID encoding */
    encodeShift = __bsr((unsigned int)maxPrimsPerGroup) + 1;
    assert( ((unsigned int)1 << encodeShift) > maxPrimsPerGroup);

    encodeMask = ((size_t)1 << encodeShift)-1;
    size_t maxGroups = ((size_t)1 << (31-encodeShift))-1;

#if 0
    DBG_PRINT(numGroups);
    DBG_PRINT(maxPrimsPerGroup);
    DBG_PRINT(numPrimitives);
    DBG_PRINT(encodeMask);
    DBG_PRINT(encodeShift);
    DBG_PRINT(maxGroups);
    DBG_PRINT(size_morton);
#endif
    if (maxPrimsPerGroup > encodeMask || numGroups > maxGroups)
    {
      DBG_PRINT(numGroups);
      DBG_PRINT(numPrimitives);
      DBG_PRINT(maxPrimsPerGroup);
      DBG_PRINT(encodeMask);
      DBG_PRINT(maxGroups);
      unsigned int primIDEncodingBits   = encodeShift;
      unsigned int groupIDEncodingBits = __bsr((unsigned int)numGroups) + 1;
      DBG_PRINT( primIDEncodingBits );
      DBG_PRINT( groupIDEncodingBits );
      FATAL("ENCODING ERROR: primIDEncodingBits + groupIDEncodingBits > 32");      
    }
  }

  void BVH4iBuilderMorton::allocateData(size_t threadCount)
  {
    /* preallocate arrays */
    const size_t additional_size = 16 * CACHELINE_SIZE;
    if (numPrimitivesOld != numPrimitives)
    {
      DBG(
	  DBG_PRINT( numPrimitivesOld );
	  DBG_PRINT( numPrimitives );
	  );

      numPrimitivesOld = numPrimitives;
      /* free previously allocated memory */
      if (morton) {
	assert(size_morton > 0);
	os_free(morton,size_morton);
      }
      if (node) { 
	assert(size_node > 0);
	os_free(node  ,size_node);
      }
      if (accel ) {
	assert(size_accel > 0);
	os_free(accel ,size_accel);
      }
      
      /* allocated memory for primrefs,nodes, and accel */
      const size_t minAllocNodes = (threadCount+1) * 2* ALLOCATOR_NODE_BLOCK_SIZE;


      const size_t numPrims      = numPrimitives+4;
      const size_t numNodes      = (size_t)((numPrimitives+3)/4);


      const size_t sizeNodeInBytes = sizeof(BVH4i::Node);
      const size_t sizeAccelInBytes = sizeof(Triangle1);


      const size_t size_morton_tmp = numPrims * sizeof(MortonID32Bit) + additional_size;

      size_node         = (double)((numNodes * MORTON_BVH4I_NODE_PREALLOC_FACTOR + minAllocNodes) * sizeNodeInBytes + additional_size) * g_memory_preallocation_factor;
      size_accel        = numPrims * sizeAccelInBytes + additional_size;
      numAllocatedNodes = size_node / sizeof(BVH4i::Node);


      DBG(DBG_PRINT(size_morton_tmp));
      DBG(DBG_PRINT(size_node));
      DBG(DBG_PRINT(size_accel));

      morton = (MortonID32Bit* ) os_malloc(size_morton_tmp); 
      node   = (BVH4i::Node*)    os_malloc(size_node  );     
      accel  = (Triangle1*)      os_malloc(size_accel );     

      memset((char*)accel + numPrims * sizeAccelInBytes,0,additional_size); // clear out as a 4-wide access is possible

      assert(morton != 0);
      assert(node   != 0);
      assert(accel  != 0);

      // memset(morton,0,size_morton_tmp);
      // memset(node  ,0,size_node);
      // memset(accel ,0,size_accel);	

      size_morton = size_morton_tmp;


      DBG(
	  DBG_PRINT( minAllocNodes );
	  DBG_PRINT( numNodes );
	  );

    }

    bvh->accel = accel;
    bvh->qbvh  = (BVH4i::Node*)node;
    bvh->size_node  = size_node;
    bvh->size_accel = size_accel;

    DBG(
	DBG_PRINT(bvh->size_node);
	DBG_PRINT(bvh->size_accel);
	DBG_PRINT(numAllocatedNodes);
	);

  }

  void BVH4iBuilderMorton::build(size_t threadIndex, size_t threadCount) 
  {    
    if (threadIndex != 0) {
      FATAL("threadIndex != 0");
    }

    if (unlikely(g_verbose >= 2))
      {
	std::cout << "building BVH4i with 32-bit Morton builder (MIC)... " << std::flush;
      }
    
    /* do some global inits first */
    initEncodingAllocateData();

    if (likely(numPrimitives == 0))
      {
	DBG(std::cout << "EMPTY SCENE BUILD" << std::endl);
	bvh->root = BVH4i::invalidNode;
	bvh->bounds = empty;
	bvh->qbvh = NULL;
	bvh->accel = NULL;
	return;
      }

    /* allocate memory arrays */
    allocateData(threadCount);

#if defined(PROFILE)
    size_t numTotalPrimitives = numPrimitives;
    std::cout << "STARTING PROFILE MODE" << std::endl << std::flush;
    std::cout << "primitives = " << numTotalPrimitives << std::endl;

    double dt_min = pos_inf;
    double dt_avg = 0.0f;
    double dt_max = neg_inf;
    size_t iterations = PROFILE_ITERATIONS;
    for (size_t i=0; i<iterations; i++) 
    {
      build_main(threadIndex,threadCount);

      dt_min = min(dt_min,dt);
      dt_avg = dt_avg + dt;
      dt_max = max(dt_max,dt);
    }
    dt_avg /= double(iterations);

    std::cout << "[DONE]" << std::endl;
    std::cout << "  min = " << 1000.0f*dt_min << "ms (" << numTotalPrimitives/dt_min*1E-6 << " Mtris/s)" << std::endl;
    std::cout << "  avg = " << 1000.0f*dt_avg << "ms (" << numTotalPrimitives/dt_avg*1E-6 << " Mtris/s)" << std::endl;
    std::cout << "  max = " << 1000.0f*dt_max << "ms (" << numTotalPrimitives/dt_max*1E-6 << " Mtris/s)" << std::endl;
    std::cout << BVH4iStatistics<BVH4i::Node>(bvh).str();

#else
    DBG(DBG_PRINT(numPrimitives));


    if (likely(numPrimitives > SINGLE_THREADED_BUILD_THRESHOLD && threadCount > 1))
      {
	DBG(std::cout << "PARALLEL BUILD" << std::endl << std::flush);
	build_main(threadIndex,threadCount);

      }
    else
      {
	/* number of primitives is small, just use single threaded mode */
	DBG(std::cout << "SERIAL BUILD" << std::endl << std::flush);
	build_main(0,1);
      }

    if (g_verbose >= 2) {
      double perf = numPrimitives/dt*1E-6;
      std::cout << "[DONE] " << 1000.0f*dt << "ms (" << perf << " Mtris/s), primitives " << numPrimitives << std::endl;
      std::cout << BVH4iStatistics<BVH4i::Node>(bvh).str();
    }
#endif
    
  }

    
  // =======================================================================================================
  // =======================================================================================================
  // =======================================================================================================

  void BVH4iBuilderMorton::initThreadState(const size_t threadID, const size_t numThreads)
  {
    const size_t numBlocks = (numPrimitives+NUM_MORTON_IDS_PER_BLOCK-1) / NUM_MORTON_IDS_PER_BLOCK;
    const size_t startID   =      ((threadID+0)*numBlocks/numThreads) * NUM_MORTON_IDS_PER_BLOCK;
    const size_t endID     = min( ((threadID+1)*numBlocks/numThreads) * NUM_MORTON_IDS_PER_BLOCK ,numPrimitives) ;
    
    assert(startID % NUM_MORTON_IDS_PER_BLOCK == 0);

    /* find first group containing startID */
    size_t group = 0, skipped = 0;
    for (; group<numGroups; group++) 
    {       
      if (unlikely(scene->get(group) == NULL)) continue;
      if (scene->get(group)->type != TRIANGLE_MESH) continue;
      const TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(group);
      if (unlikely(!mesh->isEnabled())) continue;
      if (unlikely(mesh->numTimeSteps != 1)) continue;

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
    const size_t numBlocks = (numPrimitives+NUM_MORTON_IDS_PER_BLOCK-1) / NUM_MORTON_IDS_PER_BLOCK;
    const size_t startID   =      ((threadID+0)*numBlocks/numThreads) * NUM_MORTON_IDS_PER_BLOCK;
    const size_t endID     = min( ((threadID+1)*numBlocks/numThreads) * NUM_MORTON_IDS_PER_BLOCK ,numPrimitives) ;
    assert(startID % NUM_MORTON_IDS_PER_BLOCK == 0);

    __aligned(64) Centroid_Scene_AABB bounds;
    bounds.reset();

    size_t currentID = startID;

    size_t startGroup = thread_startGroup[threadID];
    size_t offset = thread_startGroupOffset[threadID];

    mic_f bounds_centroid_min((float)pos_inf);
    mic_f bounds_centroid_max((float)neg_inf);

    for (size_t group = startGroup; group<numGroups; group++) 
    {       
      if (unlikely(scene->get(group) == NULL)) continue;
      if (unlikely(scene->get(group)->type != TRIANGLE_MESH)) continue;
      const TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(group);
      if (unlikely(!mesh->isEnabled())) continue;
      if (unlikely(mesh->numTimeSteps != 1)) continue;

      if (offset < mesh->numTriangles)
	{
	  const char* __restrict__ cptr_tri = (char*)&mesh->triangle(offset);
	  const unsigned int stride = mesh->triangles.getBufferStride();
      
	  for (size_t i=offset; i<mesh->numTriangles && currentID < endID; i++, currentID++,cptr_tri+=stride)	 
	    {
	      const TriangleMesh::Triangle& tri = *(TriangleMesh::Triangle*)cptr_tri;
	      prefetch<PFHINT_L1>(&tri + L1_PREFETCH_ITEMS);
	      prefetch<PFHINT_L2>(&tri + L2_PREFETCH_ITEMS);

	      assert( tri.v[0] < mesh->numVertices );
	      assert( tri.v[1] < mesh->numVertices );
	      assert( tri.v[2] < mesh->numVertices );

	      const mic3f v = mesh->getTriangleVertices<PFHINT_L2>(tri);

	      const mic_f bmin  = min(min(v[0],v[1]),v[2]);
	      const mic_f bmax  = max(max(v[0],v[1]),v[2]);

	      const mic_f centroid2 = bmin+bmax;
	      bounds_centroid_min = min(bounds_centroid_min,centroid2);
	      bounds_centroid_max = max(bounds_centroid_max,centroid2);
	    }
	}      
      if (unlikely(currentID == endID)) break;
      offset = 0;
    }

    store4f(&bounds.centroid2.lower,bounds_centroid_min);
    store4f(&bounds.centroid2.upper,bounds_centroid_max);
    
    global_bounds.extend_centroid_bounds_atomic(bounds); 
  }

  void BVH4iBuilderMorton::computeMortonCodes(const size_t threadID, const size_t numThreads)
  {
    const size_t numBlocks = (numPrimitives+NUM_MORTON_IDS_PER_BLOCK-1) / NUM_MORTON_IDS_PER_BLOCK;
    const size_t startID   =      ((threadID+0)*numBlocks/numThreads) * NUM_MORTON_IDS_PER_BLOCK;
    const size_t endID     = min( ((threadID+1)*numBlocks/numThreads) * NUM_MORTON_IDS_PER_BLOCK ,numPrimitives) ;
    assert(startID % NUM_MORTON_IDS_PER_BLOCK == 0);

    /* store the morton codes in 'morton' memory */
    MortonID32Bit* __restrict__ dest = ((MortonID32Bit*)morton) + startID; 

    /* compute mapping from world space into 3D grid */
    const mic_f base     = broadcast4to16f((float*)&global_bounds.centroid2.lower);
    const mic_f diagonal = 
      broadcast4to16f((float*)&global_bounds.centroid2.upper) - 
      broadcast4to16f((float*)&global_bounds.centroid2.lower);
    const mic_f scale    = select(diagonal != 0, rcp(diagonal) * mic_f(LATTICE_SIZE_PER_DIM * 0.99f),mic_f(0.0f));

    size_t currentID = startID;
    size_t offset = thread_startGroupOffset[threadID];


    mic_i mID      = mic_i::zero();
    mic_i binID3_x = mic_i::zero();
    mic_i binID3_y = mic_i::zero();
    mic_i binID3_z = mic_i::zero();

    size_t slot = 0;

    for (size_t group = thread_startGroup[threadID]; group<numGroups; group++) 
    {       
      if (unlikely(scene->get(group) == NULL)) continue;
      if (unlikely(scene->get(group)->type != TRIANGLE_MESH)) continue;
      const TriangleMesh* const mesh = scene->getTriangleMesh(group);
      if (unlikely(!mesh->isEnabled())) continue;
      if (unlikely(mesh->numTimeSteps != 1)) continue;

      const size_t numTriangles = min(mesh->numTriangles-offset,endID-currentID);
       
      const unsigned int groupCode = (group << encodeShift);

      if (offset < mesh->numTriangles)
	{
	  const char* __restrict__ cptr_tri = (char*)&mesh->triangle(offset);
	  const unsigned int stride = mesh->triangles.getBufferStride();
      
	  for (size_t i=0; i<numTriangles; i++,cptr_tri+=stride)	  
	    {
	      //const TriangleMesh::Triangle& tri = mesh->triangle(offset+i);
	      const TriangleMesh::Triangle& tri = *(TriangleMesh::Triangle*)cptr_tri;

	      prefetch<PFHINT_NT>(&tri + 16);

	      const mic3f v = mesh->getTriangleVertices<PFHINT_L2>(tri);
	      const mic_f bmin  = min(min(v[0],v[1]),v[2]);
	      const mic_f bmax  = max(max(v[0],v[1]),v[2]);

	      const mic_f cent  = bmin+bmax;
	      const mic_i binID = mic_i((cent-base)*scale);

	      mID[2*slot+1] = groupCode | (offset+i);
	      compactustore16i_low(0x1,&binID3_x[2*slot+0],binID); // extract
	      compactustore16i_low(0x2,&binID3_y[2*slot+0],binID);
	      compactustore16i_low(0x4,&binID3_z[2*slot+0],binID);
	      slot++;
	      if (unlikely(slot == NUM_MORTON_IDS_PER_BLOCK))
		{
		  const mic_i code  = bitInterleave(binID3_x,binID3_y,binID3_z);
		  const mic_i final = select(0x5555,code,mID);      
		  assert((size_t)dest % 64 == 0);
		  store16i_ngo(dest,final);	    
		  slot = 0;
		  dest += 8;
		}
	      currentID++;
	    }
	}
      offset = 0;
      if (currentID == endID) break;
    }

    if (unlikely(slot != 0))
      {
	const mic_i code  = bitInterleave(binID3_x,binID3_y,binID3_z);
	const mic_i final = select(0x5555,code,mID);      
	assert((size_t)dest % 64 == 0);
	store16i_ngo(dest,final);	    
      }
  }
  
  void BVH4iBuilderMorton::recreateMortonCodes(SmallBuildRecord& current) const
  {
    const size_t items  = current.size();
    const size_t blocks = items / NUM_MORTON_IDS_PER_BLOCK;
    const size_t rest   = items % NUM_MORTON_IDS_PER_BLOCK;

    MortonID32Bit *__restrict__ m = &morton[current.begin];

    mic_f bounds_centroid_min((float)pos_inf);
    mic_f bounds_centroid_max((float)neg_inf);


    for (size_t i=0; i<blocks; i++,m+=NUM_MORTON_IDS_PER_BLOCK)
      {
	prefetch<PFHINT_L1EX>(&morton[i+  NUM_MORTON_IDS_PER_BLOCK]);
	prefetch<PFHINT_L2EX>(&morton[i+2*NUM_MORTON_IDS_PER_BLOCK]);

	for (size_t j=0;j<NUM_MORTON_IDS_PER_BLOCK;j++)
	  {
	    const unsigned int index  = m[j].index;
	    const unsigned int primID = index & encodeMask; 
	    const unsigned int geomID = index >> encodeShift; 
	    const TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
	    const TriangleMesh::Triangle& tri = mesh->triangle(primID);

	    const mic3f v = mesh->getTriangleVertices<PFHINT_L1>(tri);
	    const mic_f bmin  = min(min(v[0],v[1]),v[2]);
	    const mic_f bmax  = max(max(v[0],v[1]),v[2]);

	    const mic_f centroid2 = bmin+bmax;
	    bounds_centroid_min = min(bounds_centroid_min,centroid2);
	    bounds_centroid_max = max(bounds_centroid_max,centroid2);
	  }
      }   

    if (rest)
      {
	for (size_t j=0;j<rest;j++)
	  {
	    const unsigned int index  = m[j].index;
	    const unsigned int primID = index & encodeMask; 
	    const unsigned int geomID = index >> encodeShift; 

	    const TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
	    const TriangleMesh::Triangle& tri = mesh->triangle(primID);

	    const mic3f v = mesh->getTriangleVertices<PFHINT_L1>(tri);
	    const mic_f bmin  = min(min(v[0],v[1]),v[2]);
	    const mic_f bmax  = max(max(v[0],v[1]),v[2]);

	    const mic_f centroid2 = bmin+bmax;
	    bounds_centroid_min = min(bounds_centroid_min,centroid2);
	    bounds_centroid_max = max(bounds_centroid_max,centroid2);
	  }
      }


    const mic_f base     = bounds_centroid_min;
    const mic_f diagonal = bounds_centroid_max - bounds_centroid_min;
    const mic_f scale    = select(diagonal != 0,rcp(diagonal) * mic_f(LATTICE_SIZE_PER_DIM * 0.99f),mic_f(0.0f));
    
    mic_i binID3_x = mic_i::zero();
    mic_i binID3_y = mic_i::zero();
    mic_i binID3_z = mic_i::zero();

    m = &morton[current.begin];

    for (size_t i=0; i<blocks; i++,m+=NUM_MORTON_IDS_PER_BLOCK)
      {
	for (size_t j=0;j<NUM_MORTON_IDS_PER_BLOCK;j++)
	  {
	    const unsigned int index  = m[j].index;
	    const unsigned int primID = index & encodeMask; 
	    const unsigned int geomID = index >> encodeShift; 

	    const TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
	    const TriangleMesh::Triangle& tri = mesh->triangle(primID);

	    const mic3f v = mesh->getTriangleVertices(tri);
	    const mic_f bmin  = min(min(v[0],v[1]),v[2]);
	    const mic_f bmax  = max(max(v[0],v[1]),v[2]);
	    const mic_f centroid = bmin+bmax;
	    const mic_i binID = mic_i((centroid-base)*scale);

	    compactustore16i_low(0x1,&binID3_x[2*j+0],binID); 
	    compactustore16i_low(0x2,&binID3_y[2*j+0],binID);
	    compactustore16i_low(0x4,&binID3_z[2*j+0],binID);	    
	  }

	const mic_i mID = uload16i((int*)m);
	const mic_i code  = bitInterleave(binID3_x,binID3_y,binID3_z);
	const mic_i final = select(0x5555,code,mID);      
	ustore16i(m,final);	
      }
    if (rest)
      {
	for (size_t j=0;j<rest;j++)
	  {
	    const unsigned int index  = m[j].index;
	    const unsigned int primID = index & encodeMask; 
	    const unsigned int geomID = index >> encodeShift; 

	    const TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
	    const TriangleMesh::Triangle& tri = mesh->triangle(primID);

	    const mic3f v = mesh->getTriangleVertices(tri);
	    const mic_f bmin  = min(min(v[0],v[1]),v[2]);
	    const mic_f bmax  = max(max(v[0],v[1]),v[2]);
	    const mic_f centroid = bmin+bmax;
	    const mic_i binID = mic_i((centroid-base)*scale);

	    compactustore16i_low(0x1,&binID3_x[2*j+0],binID); 
	    compactustore16i_low(0x2,&binID3_y[2*j+0],binID);
	    compactustore16i_low(0x4,&binID3_z[2*j+0],binID);	    
	  }
	const mic_m mask = ((unsigned int)1 << (2*rest))-1;
	const mic_i mID = uload16i((int*)m);
	const mic_i code  = bitInterleave(binID3_x,binID3_y,binID3_z);
	const mic_i final = select(0x5555,code,mID);      
	compactustore16i(mask,(int*)m,final);		
      }       

    quicksort_insertionsort_ascending<MortonID32Bit,32>(morton,current.begin,current.end-1); 

#if defined(DEBUG)
    for (size_t i=current.begin; i<current.end-1; i++)
      assert(morton[i].code <= morton[i+1].code);
#endif	    

  }


  void BVH4iBuilderMorton::radixsort(const size_t threadID, const size_t numThreads)
  {
    const size_t numBlocks = (numPrimitives+NUM_MORTON_IDS_PER_BLOCK-1) / NUM_MORTON_IDS_PER_BLOCK;
    const size_t startID   = ((threadID+0)*numBlocks/numThreads) * NUM_MORTON_IDS_PER_BLOCK;
    const size_t endID     = ((threadID+1)*numBlocks/numThreads) * NUM_MORTON_IDS_PER_BLOCK;
    assert(startID % NUM_MORTON_IDS_PER_BLOCK == 0);
    assert(endID % NUM_MORTON_IDS_PER_BLOCK == 0);

    assert(((numThreads)*numBlocks/numThreads) * NUM_MORTON_IDS_PER_BLOCK == ((numPrimitives+7)&(-8)));

    MortonID32Bit* __restrict__ mortonID[2];
    mortonID[0] = (MortonID32Bit*) morton; 
    mortonID[1] = (MortonID32Bit*) node;


    /* we need 4 iterations to process all 32 bits */
    for (size_t b=0; b<4; b++)
    {
      const MortonID32Bit* __restrict__ const src = (MortonID32Bit*)mortonID[((b+0)%2)];
      MortonID32Bit*       __restrict__ const dst = (MortonID32Bit*)mortonID[((b+1)%2)];

      __assume_aligned(&radixCount[threadID][0],64);
      
      /* count how many items go into the buckets */

#pragma unroll(16)
      for (size_t i=0; i<16; i++)
	store16i(&radixCount[threadID][i*16],mic_i::zero());


      for (size_t i=startID; i<endID; i+=NUM_MORTON_IDS_PER_BLOCK) {
	prefetch<PFHINT_NT>(&src[i+L1_PREFETCH_ITEMS]);
	prefetch<PFHINT_L2>(&src[i+L2_PREFETCH_ITEMS]);
	
#pragma unroll(NUM_MORTON_IDS_PER_BLOCK)
	for (unsigned long j=0;j<NUM_MORTON_IDS_PER_BLOCK;j++)
	  {
	    const unsigned int index = src[i+j].getByte(b);
	    radixCount[threadID][index]++;
	  }
      }

      scene->lockstep_scheduler.syncThreads(threadID,numThreads);


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
      for (size_t i=startID; i<endID; i+=NUM_MORTON_IDS_PER_BLOCK) {
	prefetch<PFHINT_NT>(&src[i+L1_PREFETCH_ITEMS]);
	prefetch<PFHINT_L2>(&src[i+L2_PREFETCH_ITEMS]);

#pragma nounroll
	for (unsigned long j=0;j<NUM_MORTON_IDS_PER_BLOCK;j++)
	  {
	    const unsigned int index = src[i+j].getByte(b);
	    assert(index < RADIX_BUCKETS);
	    dst[offset[index]] = src[i+j];
	    prefetch<PFHINT_L2EX>(&dst[offset[index]+L1_PREFETCH_ITEMS]);
	    offset[index]++;
	  }
	evictL2(&src[i]);
      }

      if (b<3) scene->lockstep_scheduler.syncThreads(threadID,numThreads);

    }
  }

  void BVH4iBuilderMorton::createTopLevelTree(const size_t threadID, const size_t numThreads)
  {
    size_t taskID = threadID;
    __aligned(64) SmallBuildRecord children[BVH4i::N];

    
    while(taskID < numBuildRecords)
      {
	SmallBuildRecord &sbr = buildRecords[taskID];
	if (sbr.size() > topLevelItemThreshold)
	  {
	    const size_t numChildren = createSingleBVH4iNode(sbr,children);
	    buildRecords[taskID] = children[0];
	    if (numChildren > 1)
	      {
		const unsigned int dest = numBuildRecordCounter.add(numChildren-1);
		for (size_t i=0;i<numChildren-1;i++)
		  buildRecords[dest + i] = children[i+1];
	      }
	  }
	taskID += numThreads;
      }   
  }

  void BVH4iBuilderMorton::recurseSubMortonTrees(const size_t threadID, const size_t numThreads)
  {
    NodeAllocator alloc(atomicID,numAllocatedNodes);
    
    while (true)
    {
      const unsigned int taskID = scene->lockstep_scheduler.taskCounter.inc();
      if (taskID >= numBuildRecords) break;
      
      SmallBuildRecord &br = buildRecords[taskID];
      
      BBox3fa bounds = recurse(br,alloc);     
      
      /* mark toplevel of tree */
      node[br.parentNodeID].setBounds(br.parentLocalID,bounds);
      node[br.parentNodeID].upper[br.parentLocalID].child = BVH4I_TOP_LEVEL_MARKER;
    }    

  }
  
  // =======================================================================================================
  // =======================================================================================================
  // =======================================================================================================


  void BVH4iBuilderMorton::split_fallback(SmallBuildRecord& current, SmallBuildRecord& leftChild, SmallBuildRecord& rightChild)
  {
    unsigned int blocks4 = (current.items()+3)/4;
    unsigned int center = current.begin + (blocks4/2)*4; 

    assert(center != current.begin);
    assert(center != current.end);

    leftChild.init(current.begin,center);
    rightChild.init(center,current.end);
  }
		

  __forceinline BBox3fa BVH4iBuilderMorton::createSmallLeaf(SmallBuildRecord& current)
  {    
    assert(current.size() > 0);
    mic_f bounds_min(pos_inf);
    mic_f bounds_max(neg_inf);

    Vec3fa lower(pos_inf);
    Vec3fa upper(neg_inf);
    size_t items = current.size();
    size_t start = current.begin;
    assert(items<=4);

    const mic_i morton_mask(encodeMask);
    const mic_i morton_shift(encodeShift);

    prefetch<PFHINT_L2>(&morton[start+8]);

    for (size_t i=0; i<items; i++) 
      {	
	const unsigned int index = morton[start+i].index;
	const unsigned int primID = index & encodeMask; 
	const unsigned int geomID = index >> encodeShift; 


	const mic_i morton_index(morton[start+i].index);
	const mic_i morton_primID = morton_index & morton_mask;
	const mic_i morton_geomID = morton_index >> morton_shift;

	const TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
	const TriangleMesh::Triangle& tri = mesh->triangle(primID);

	const mic3f v = mesh->getTriangleVertices(tri);
	const mic_f v0 = v[0];
	const mic_f v1 = v[1];
	const mic_f v2 = v[2];

	//WARNING: zero last component

	const mic_f tri_accel = initTriangle1(v0,v1,v2,morton_geomID,morton_primID,mic_i(mesh->mask));

	bounds_min = min(bounds_min,min(v0,min(v1,v2)));
	bounds_max = max(bounds_max,max(v0,max(v1,v2)));
	store16f_ngo(&accel[start+i],tri_accel);
      }

    store3f(&node[current.parentNodeID].lower[current.parentLocalID],bounds_min);
    store3f(&node[current.parentNodeID].upper[current.parentLocalID],bounds_max);
    createBVH4iLeaf(node[current.parentNodeID].lower[current.parentLocalID].child,start,items);
    __aligned(64) BBox3fa bounds;
    store4f(&bounds.lower,bounds_min);
    store4f(&bounds.upper,bounds_max);

    return bounds;
  }


  BBox3fa BVH4iBuilderMorton::createLeaf(SmallBuildRecord& current, NodeAllocator& alloc)
  {
#if defined(DEBUG)
    if (current.depth > BVH4i::maxBuildDepthLeaf) 
      THROW_RUNTIME_ERROR("ERROR: depth limit reached");
#endif
    
    /* create leaf for few primitives */
    if (current.size() <= MORTON_LEAF_THRESHOLD) {     
      return createSmallLeaf(current);
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
    const size_t currentIndex = alloc.get(1);
   
    /* init used/unused nodes */
    mic_f init_lower = broadcast4to16f(&BVH4i::initQBVHNode[0]);
    mic_f init_upper = broadcast4to16f(&BVH4i::initQBVHNode[1]);

    store16f_ngo((float*)&node[currentIndex].lower,init_lower);
    store16f_ngo((float*)&node[currentIndex].upper,init_upper);


    __aligned(64) BBox3fa bounds; 
    bounds = empty;
    /* recurse into each child */
    for (size_t i=0; i<numChildren; i++) {
      children[i].parentNodeID  = currentIndex;
      children[i].parentLocalID = i;
      children[i].depth = current.depth+1;
      bounds.extend( createLeaf(children[i],alloc) );
    }

    store3f(&node[current.parentNodeID].lower[current.parentLocalID],broadcast4to16f(&bounds.lower));
    store3f(&node[current.parentNodeID].upper[current.parentLocalID],broadcast4to16f(&bounds.upper));

    createBVH4iNode<4>(node[current.parentNodeID].lower[current.parentLocalID].child,currentIndex);

    return bounds;
  }  

  __forceinline bool BVH4iBuilderMorton::split(SmallBuildRecord& current,
                                               SmallBuildRecord& left,
                                               SmallBuildRecord& right)
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
    unsigned int begin = current.begin;
    unsigned int end   = current.end;
    while (begin + 1 != end) {
      const unsigned int mid = (begin+end)/2;
      const unsigned int bit = morton[mid].code & bitmask;
      if (bit == 0) begin = mid; else end = mid;
    }
    unsigned int center = end;
#if defined(DEBUG)      
    for (unsigned int i=begin;  i<center; i++) assert((morton[i].code & bitmask) == 0);
    for (unsigned int i=center; i<end;    i++) assert((morton[i].code & bitmask) == bitmask);
#endif
    
    left.init(current.begin,center);
    right.init(center,current.end);
    return true;
  }

  size_t BVH4iBuilderMorton::createSingleBVH4iNode(SmallBuildRecord& current, SmallBuildRecord *__restrict__ const children)
  {

    /* create leaf node */
    if (unlikely(current.size() <= BVH4iBuilderMorton::MORTON_LEAF_THRESHOLD)) {
      children[0] = current;
      return 1;
    }

    /* fill all 4 children by always splitting the one with the largest number of primitives */
    __assume_aligned(children,sizeof(SmallBuildRecord));

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
      children[0] = current;
      return 1;
    }

    /* allocate next four nodes and prefetch them */
    const size_t currentIndex = allocGlobalNode(1);    

    mic_f init_lower = broadcast4to16f(&BVH4i::initQBVHNode[0]);
    mic_f init_upper = broadcast4to16f(&BVH4i::initQBVHNode[1]);

    store16f_ngo((float*)&node[currentIndex].lower,init_lower);
    store16f_ngo((float*)&node[currentIndex].upper,init_upper);

    /* recurse into each child */
    for (size_t i=0; i<numChildren; i++) 
      {
	children[i].parentNodeID  = currentIndex;
	children[i].parentLocalID = i;
      }

    createBVH4iNode<4>(node[current.parentNodeID].lower[current.parentLocalID].child,currentIndex);

    return numChildren;
  }

  
  BBox3fa BVH4iBuilderMorton::recurse(SmallBuildRecord& current, 
				     NodeAllocator& alloc) 
  {
    assert(current.size() > 0);

    __aligned(64) SmallBuildRecord children[BVH4i::N];

    /* create leaf node */
    if (unlikely(current.size() <= BVH4iBuilderMorton::MORTON_LEAF_THRESHOLD)) {
      return createSmallLeaf(current);
    }
    if (unlikely(current.depth >= BVH4i::maxBuildDepth)) {
      return createLeaf(current,alloc); 
    }

    /* fill all 4 children by always splitting the one with the largest number of primitives */
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
      return createSmallLeaf(current);
    }

    /* allocate next four nodes and prefetch them */
    const size_t currentIndex = alloc.get(1);    

    /* init used/unused nodes */
    mic_f init_lower = broadcast4to16f(&BVH4i::initQBVHNode[0]);
    mic_f init_upper = broadcast4to16f(&BVH4i::initQBVHNode[1]);

    store16f_ngo((float*)&node[currentIndex].lower,init_lower);
    store16f_ngo((float*)&node[currentIndex].upper,init_upper);

    /* recurse into each child */
    __aligned(64) BBox3fa bounds;
    bounds = empty;
    size_t child_nodes = 0;
    for (size_t i=0; i<numChildren; i++) 
    {
      children[i].parentNodeID = currentIndex;
      children[i].parentLocalID = i;

      if (children[i].size() <= BVH4iBuilderMorton::MORTON_LEAF_THRESHOLD)
	{
	  bounds.extend( createSmallLeaf(children[i]) );
	}
      else
	{
	  bounds.extend( recurse(children[i],alloc) ); // recurse 
	}
    }

    store3f(&node[current.parentNodeID].lower[current.parentLocalID],broadcast4to16f(&bounds.lower));
    store3f(&node[current.parentNodeID].upper[current.parentLocalID],broadcast4to16f(&bounds.upper));
    createBVH4iNode<4>(node[current.parentNodeID].lower[current.parentLocalID].child,currentIndex);

    if (enableTreeRotations)
      {
	BVH4i::NodeRef local_root = node[current.parentNodeID].child(current.parentLocalID);
	BVH4iRotate::rotateNode(bvh,local_root);
      }

    return bounds;
  }

  BBox3fa BVH4iBuilderMorton::refit(const BVH4i::NodeRef &ref)
  {    
    if (unlikely(ref.isLeaf()))
      {
	return BBox3fa( empty );
      }

    BVH4i::Node *n = (BVH4i::Node*)ref.node(node);

    BBox3fa parentBounds = empty;

    for (size_t i=0;i<BVH4i::N;i++)
      {
	if (n->child(i) == BVH4i::invalidNode) break;
	
	if (n->child(i).isLeaf())
	  {
	    parentBounds.extend( n->bounds(i) );
	  }
	else
	  {
	    BBox3fa bounds = refit( n->child(i) );

	    n->setBounds(i,bounds);
	    parentBounds.extend( bounds );
	  }
      }
    return parentBounds;
  }    

  BBox3fa BVH4iBuilderMorton::refit_toplevel(const BVH4i::NodeRef &ref)
  {    
    if (unlikely(ref.isLeaf()))
      {
	return BBox3fa( empty );
      }

    BVH4i::Node *n = (BVH4i::Node*)ref.node(node);

    BBox3fa parentBounds = empty;

    for (size_t i=0;i<BVH4i::N;i++)
      {
	if (n->child(i) == BVH4i::invalidNode) break;
	
	if (n->child(i).isLeaf())
	  parentBounds.extend( n->bounds(i) );
	else if (n->upper[i].child == BVH4I_TOP_LEVEL_MARKER)
	  {
	    BVH4i::Node *c = (BVH4i::Node*)n->child(i).node(node);
	    parentBounds.extend( c->bounds() );	    
	    n->setBounds(i,c->bounds() );
	  }
	else
	  {
	    BBox3fa bounds = refit_toplevel( n->child(i) );
	    n->setBounds(i,bounds);
	    parentBounds.extend( bounds );
	  }
      }

    if (enableTreeRotations)
      BVH4iRotate::rotateNode(bvh,ref);

    return parentBounds;

  }

  void BVH4iBuilderMorton::build_main (const size_t threadIndex, const size_t threadCount)
  { 
    DBG(PING);

    /* start measurement */
    double t0 = 0.0f;

#if !defined(PROFILE)
    if (g_verbose >= 2) 
#endif
      t0 = getSeconds();

    TIMER(std::cout << std::endl);
    TIMER(double msec = 0.0);

    /* init thread state */
    TIMER(msec = getSeconds());
    scene->lockstep_scheduler.dispatchTask( task_initThreadState, this, threadIndex, threadCount );
    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "task_initThreadState " << 1000. * msec << " ms" << std::endl << std::flush);

    /* compute scene bounds */
    TIMER(msec = getSeconds());
    global_bounds.reset();
    scene->lockstep_scheduler.dispatchTask( task_computeBounds, this, threadIndex, threadCount );
    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "task_computeBounds " << 1000. * msec << " ms" << std::endl << std::flush);
    TIMER(DBG_PRINT(global_bounds));
    
    /* compute morton codes */
    TIMER(msec = getSeconds());
    scene->lockstep_scheduler.dispatchTask( task_computeMortonCodes, this, threadIndex, threadCount );   

    /* padding */
    MortonID32Bit* __restrict__ const dest = (MortonID32Bit*)morton;
    
    for (size_t i=numPrimitives; i<((numPrimitives+7)&(-8)); i++) {
      dest[i].code  = 0xffffffff; 
      dest[i].index = 0;
    }

    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "task_computeMortonCodes " << 1000. * msec << " ms" << std::endl << std::flush);

 

    /* sort morton codes */
    TIMER(msec = getSeconds());
    scene->lockstep_scheduler.dispatchTask( task_radixsort, this, threadIndex, threadCount );

#if defined(DEBUG)
    for (size_t i=1; i<((numPrimitives+7)&(-8)); i++)
      {
	if (morton[i-1].code > morton[i].code)
	  {
	    DBG_PRINT( i );
	    DBG_PRINT( morton[i-1].code );
	    DBG_PRINT( morton[i].code );
	  }

	assert(morton[i-1].code <= morton[i].code);
      }

    for (size_t i=numPrimitives; i<((numPrimitives+7)&(-8)); i++) {
      assert(dest[i].code  == 0xffffffff); 
      assert(dest[i].index == 0);
    }

#endif	    

    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "task_radixsort " << 1000. * msec << " ms -> " << numPrimitives / msec / 1E+6 << "M key/value pairs per sec" <<  std::endl << std::flush);

    TIMER(msec = getSeconds());

    /* build and extract top-level tree */
    numBuildRecords = 0;
    atomicID.reset(BVH4i::N);
    topLevelItemThreshold = max((numPrimitives + threadCount-1)/((threadCount)),(size_t)64);

    SmallBuildRecord br;
    br.init(0,numPrimitives);
    br.parentNodeID = 0;
    br.parentLocalID = 0;
    br.depth = 1;


    buildRecords[0] = br;
    numBuildRecords = 1;
    size_t iterations = 0;
    while(numBuildRecords < threadCount*3)
      {
	numBuildRecordCounter.reset(numBuildRecords);
	scene->lockstep_scheduler.dispatchTask( task_createTopLevelTree, this, threadIndex, threadCount );
	iterations++;

	if (unlikely(numBuildRecords == numBuildRecordCounter)) { break; }

	numBuildRecords = numBuildRecordCounter;
      }

    /* sort all subtasks by size */

    quicksort_insertionsort_decending<SmallBuildRecord,16>(buildRecords,0,numBuildRecords-1);



    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "create top level " << 1000. * msec << " ms" << std::endl << std::flush);
    TIMER(DBG_PRINT(numBuildRecords));

    TIMER(msec = getSeconds());
    
    /* build sub-trees */
    scene->lockstep_scheduler.dispatchTask( task_recurseSubMortonTrees, this, threadIndex, threadCount );

    DBG(DBG_PRINT(atomicID));

    numNodes = atomicID;

    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "task_recurseSubMortonTrees " << 1000. * msec << " ms" << std::endl << std::flush);

    TIMER(msec = getSeconds());

    /* refit toplevel part of tree */
    BBox3fa rootBounds = refit_toplevel(node->child(0));

    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "refit top level " << 1000. * msec << " ms" << std::endl << std::flush);

    bvh->root   = node->child(0); 
    bvh->bounds = node->child(0).isLeaf() ? node->bounds(0) : rootBounds;

    /* stop measurement */
#if !defined(PROFILE)
    if (g_verbose >= 2) 
#endif
      dt = getSeconds()-t0;

  }

}


