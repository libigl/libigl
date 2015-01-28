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
#include "bvh4i/bvh4i_statistics.h"


#define THRESHOLD_FOR_SUBTREE_RECURSION         64
#define BUILD_RECORD_PARALLEL_SPLIT_THRESHOLD 1024
#define SINGLE_THREADED_BUILD_THRESHOLD        512
#define ENABLE_TASK_STEALING
#define ENABLE_FILL_PER_CORE_WORK_QUEUES


#define L1_PREFETCH_ITEMS 2
#define L2_PREFETCH_ITEMS 16

#define TIMER(x) 
#define DBG(x) 

//#define PROFILE
#define PROFILE_ITERATIONS 100

#define MEASURE_MEMORY_ALLOCATION_TIME 0

// TODO: CHECK     const float voxelArea    = current.cs_AABB.sceneArea();
//                 const float centroidArea = current.cs_AABB.centroidArea();

namespace embree
{
#if defined(DEBUG)
  static AtomicMutex mtx;
#endif

  static double dt = 0.0f;

  Builder* BVH4iBuilder::create (void* accel, BuildSource* source, void* geometry, size_t mode ) 
  { 
    Builder* builder = NULL;

    switch( mode )
      {

      case BVH4I_BUILDER_DEFAULT:
	builder = new BVH4iBuilder((BVH4i*)accel,source,geometry);
	break;

      case BVH4I_BUILDER_PRESPLITS:
	builder = new BVH4iBuilderPreSplits((BVH4i*)accel,source,geometry);
	break;

      case BVH4I_BUILDER_VIRTUAL_GEOMETRY:
	builder = new BVH4iBuilderVirtualGeometry((BVH4i*)accel,source,geometry);
	break;

      default:
	throw std::runtime_error("ERROR: unknown BVH4iBuilder mode selected");	
      }
    return builder;
  }

  // =======================================================================================================
  // =======================================================================================================
  // =======================================================================================================


  BVH4iBuilder::BVH4iBuilder (BVH4i* bvh, BuildSource* source, void* geometry)
    : source(source), 
      scene((Scene*)geometry), 
      bvh(bvh), 
      numPrimitives((size_t)-1), 
      numNodes(0), 
      numAllocatedNodes(0), 
      prims(NULL), 
      node(NULL), 
      accel(NULL), 
      size_prims(0),
      numNodesToAllocate(BVH4i::N)
  {
    DBG(PING);
  }

  BVH4iBuilder::~BVH4iBuilder()
  {
    DBG(PING);
    if (prims)  {
      assert(size_prims > 0);
      os_free(prims,size_prims);
    }    
  }

  size_t BVH4iBuilder::getNumPrimitives()
  {
    /* count total number of triangles */
    size_t primitives = 0;       
    for (size_t i=0;i<scene->size();i++)
      {
	if (unlikely(scene->get(i) == NULL)) continue;
	if (unlikely((scene->get(i)->type != TRIANGLE_MESH))) continue;
	if (unlikely(!scene->get(i)->isEnabled())) continue;
	const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(i);
	if (unlikely(mesh->numTimeSteps != 1)) continue;
	primitives += mesh->numTriangles;
      }
    return primitives;	
  
  }

  void BVH4iBuilder::allocateMemoryPools(const size_t numPrims, 
					 const size_t numNodes,
					 const size_t sizeNodeInBytes,
					 const size_t sizeAccelInBytes)
  {
#if MEASURE_MEMORY_ALLOCATION_TIME == 1
    double msec = 0.0;
    msec = getSeconds();
#endif

    const size_t additional_size = 16 * CACHELINE_SIZE;

    /* free previously allocated memory */

    if (prims)  {
      assert(size_prims > 0);
      os_free(prims,size_prims);
    }
    if (node  ) {
      assert(bvh->size_node > 0);
      os_free(node ,bvh->size_node);
    }
    if (accel ) {
      assert(bvh->size_accel > 0);
      os_free(accel,bvh->size_accel);
    }
      
    // === allocated memory for primrefs,nodes, and accel ===
    const size_t size_primrefs = numPrims * sizeof(PrimRef) + additional_size;
    const size_t size_node     = numNodes * BVH_NODE_PREALLOC_FACTOR * sizeNodeInBytes + additional_size;
    const size_t size_accel    = numPrims * sizeAccelInBytes + additional_size;

    numAllocatedNodes = size_node / sizeof(BVHNode);
      
    DBG(DBG_PRINT(numAllocatedNodes));
    DBG(DBG_PRINT(size_primrefs));
    DBG(DBG_PRINT(size_node));
    DBG(DBG_PRINT(size_accel));

    prims = (PrimRef  *) os_malloc(size_primrefs); 
    node  = (BVHNode  *) os_malloc(size_node);
    accel = (Triangle1*) os_malloc(size_accel);

    assert(prims  != 0);
    assert(node   != 0);
    assert(accel  != 0);

#if 0
    memset(prims,0,size_primrefs);
    memset(node,0,size_node);
    memset(accel,0,size_accel);
#endif

    bvh->accel = accel;
    bvh->qbvh  = (BVH4i::Node*)node;
    bvh->size_node  = size_node;
    bvh->size_accel = size_accel;

    size_prims = size_primrefs;    

#if MEASURE_MEMORY_ALLOCATION_TIME == 1
    msec = getSeconds()-msec;    
    std::cout << "allocation time " << 1000. * msec << " ms" << std::endl << std::flush;
#endif
  }
  
  void BVH4iBuilder::allocateData(const size_t threadCount, const size_t totalNumPrimitives)
  {
    DBG(PING);
    size_t numPrimitivesOld = numPrimitives;
    numPrimitives = totalNumPrimitives;
    DBG(DBG_PRINT(numPrimitives));


    if (numPrimitivesOld != numPrimitives)
      {
	const size_t numPrims = numPrimitives+4;
	const size_t minAllocNodes = numPrims ? threadCount * ALLOCATOR_NODE_BLOCK_SIZE * 4: 16;
	const size_t numNodes = max((size_t)(numPrims * BVH_NODE_PREALLOC_FACTOR),minAllocNodes);
	allocateMemoryPools(numPrims,numNodes);
      }
  }

  void BVH4iBuilder::printBuilderName()
  {
    std::cout << "building BVH4i with SAH builder (MIC) ... " << std::endl;    
  }


  void BVH4iBuilder::build(size_t threadIndex, size_t threadCount) 
  {
    DBG(PING);
    const size_t totalNumPrimitives = getNumPrimitives();


    DBG(DBG_PRINT(totalNumPrimitives));

    /* print builder name */
    if (unlikely(g_verbose >= 1)) printBuilderName();

    /* allocate BVH data */
    allocateData(TaskScheduler::getNumThreads(),totalNumPrimitives);

    LockStepTaskScheduler::init(TaskScheduler::getNumThreads()); 

    if (likely(numPrimitives > SINGLE_THREADED_BUILD_THRESHOLD && TaskScheduler::getNumThreads() > 1) )
      {
	DBG(std::cout << "PARALLEL BUILD" << std::endl);

#if defined(PROFILE)

	std::cout << "STARTING PROFILE MODE" << std::endl << std::flush;
	std::cout << "primitives = " << totalNumPrimitives << std::endl;
	double dt_min = pos_inf;
	double dt_avg = 0.0f;
	double dt_max = neg_inf;
	size_t iterations = PROFILE_ITERATIONS;
	for (size_t i=0; i<iterations; i++) 
	  {
	    TaskScheduler::executeTask(threadIndex,threadCount,_build_parallel,this,TaskScheduler::getNumThreads(),"build_parallel");
	    dt_min = min(dt_min,dt);
	    dt_avg = dt_avg + dt;
	    dt_max = max(dt_max,dt);
	  }
	dt_avg /= double(iterations);

	std::cout << "[DONE]" << std::endl;
	std::cout << "  min = " << 1000.0f*dt_min << "ms (" << totalNumPrimitives/dt_min*1E-6 << " Mtris/s)" << std::endl;
	std::cout << "  avg = " << 1000.0f*dt_avg << "ms (" << totalNumPrimitives/dt_avg*1E-6 << " Mtris/s)" << std::endl;
	std::cout << "  max = " << 1000.0f*dt_max << "ms (" << totalNumPrimitives/dt_max*1E-6 << " Mtris/s)" << std::endl;
	std::cout << BVH4iStatistics(bvh).str();

#else

	TaskScheduler::executeTask(threadIndex,threadCount,_build_parallel,this,TaskScheduler::getNumThreads(),"build_parallel");
#endif
      }
    else
      {
	/* number of primitives is small, just use single threaded mode */
	if (likely(numPrimitives > 0))
	  {
	    DBG(std::cout << "SERIAL BUILD" << std::endl);
	    build_parallel(0,1,0,0,NULL);
	  }
	else
	  {
	    DBG(std::cout << "EMPTY SCENE BUILD" << std::endl);
	    /* handle empty scene */
	    for (size_t i=0;i<4;i++)
	      bvh->qbvh[0].setInvalid(i);
	    for (size_t i=0;i<4;i++)
	      bvh->qbvh[1].setInvalid(i);
	    bvh->qbvh[0].lower[0].child = BVH4i::NodeRef(128);
	    bvh->root = bvh->qbvh[0].lower[0].child; 
	    bvh->bounds = empty;
	  }
      }



    if (g_verbose >= 2) {
      double perf = totalNumPrimitives/dt*1E-6;
      std::cout << "[DONE] " << 1000.0f*dt << "ms (" << perf << " Mtris/s), primitives " << numPrimitives << std::endl;
      std::cout << BVH4iStatistics(bvh).str();
    }

  }


  // =======================================================================================================
  // =======================================================================================================
  // =======================================================================================================

  void BVH4iBuilder::computePrimRefsTriangles(const size_t threadID, const size_t numThreads) 
  {
    DBG(PING);
    const size_t numGroups = scene->size();
    const size_t startID = (threadID+0)*numPrimitives/numThreads;
    const size_t endID   = (threadID+1)*numPrimitives/numThreads;
    
    PrimRef *__restrict__ const prims     = this->prims;

    // === find first group containing startID ===
    unsigned int g=0, numSkipped = 0;
    for (; g<numGroups; g++) {       
      if (unlikely(scene->get(g) == NULL)) continue;
      if (unlikely(scene->get(g)->type != TRIANGLE_MESH)) continue;
      const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(g);
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

    unsigned int num = 0;
    unsigned int currentID = startID;
    unsigned int offset = startID - numSkipped;

    __aligned(64) PrimRef local_prims[2];
    size_t numLocalPrims = 0;
    PrimRef *__restrict__ dest = &prims[currentID];

    for (; g<numGroups; g++) 
    {
      if (unlikely(scene->get(g) == NULL)) continue;
      if (unlikely(scene->get(g)->type != TRIANGLE_MESH)) continue;
      const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(g);
      if (unlikely(!mesh->isEnabled())) continue;
      if (unlikely(mesh->numTimeSteps != 1)) continue;

      for (unsigned int i=offset; i<mesh->numTriangles && currentID < endID; i++, currentID++)	 
      { 			    
	//DBG_PRINT(currentID);
	const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(i);
	prefetch<PFHINT_L2>(&tri + L2_PREFETCH_ITEMS);
	prefetch<PFHINT_L1>(&tri + L1_PREFETCH_ITEMS);

	const float *__restrict__ const vptr0 = (float*)&mesh->vertex(tri.v[0]);
	const float *__restrict__ const vptr1 = (float*)&mesh->vertex(tri.v[1]);
	const float *__restrict__ const vptr2 = (float*)&mesh->vertex(tri.v[2]);

	const mic_f v0 = broadcast4to16f(vptr0);
	const mic_f v1 = broadcast4to16f(vptr1);
	const mic_f v2 = broadcast4to16f(vptr2);

	const mic_f bmin = min(min(v0,v1),v2);
	const mic_f bmax = max(max(v0,v1),v2);
	bounds_scene_min = min(bounds_scene_min,bmin);
	bounds_scene_max = max(bounds_scene_max,bmax);
	const mic_f centroid2 = bmin+bmax;
	bounds_centroid_min = min(bounds_centroid_min,centroid2);
	bounds_centroid_max = max(bounds_centroid_max,centroid2);

	store4f(&local_prims[numLocalPrims].lower,bmin);
	store4f(&local_prims[numLocalPrims].upper,bmax);	
	local_prims[numLocalPrims].lower.a = g;
	local_prims[numLocalPrims].upper.a = i;

	//DBG_PRINT( local_prims[numLocalPrims] );

	numLocalPrims++;
	if (unlikely(((size_t)dest % 64) != 0) && numLocalPrims == 1)
	  {
	    *dest = local_prims[0];
	    dest++;
	    numLocalPrims--;
	  }
	else
	  {
	    const mic_f twoAABBs = load16f(local_prims);
	    if (numLocalPrims == 2)
	      {
		numLocalPrims = 0;
		store16f_ngo(dest,twoAABBs);
		dest+=2;
	      }
	  }	
      }
      if (currentID == endID) break;
      offset = 0;
    }

    /* is there anything left in the local queue? */
    if (numLocalPrims % 2 != 0)
      *dest = local_prims[0];

    /* update global bounds */
    Centroid_Scene_AABB bounds;
    
    store4f(&bounds.centroid2.lower,bounds_centroid_min);
    store4f(&bounds.centroid2.upper,bounds_centroid_max);
    store4f(&bounds.geometry.lower,bounds_scene_min);
    store4f(&bounds.geometry.upper,bounds_scene_max);

    global_bounds.extend_atomic(bounds);    
  }

  __forceinline void reorderBVHNodesOnArea(BVHNode *__restrict__ const bptr)
  {
    float node_area[4];
    size_t valid = 0;
    for (size_t i=0;i<4;i++) 
      {
	node_area[i] = area( bptr[i] );
	if ( node_area[i] > 0.0f) valid++;
      }
    
    if (valid == 0) return;

    assert( valid >= 2 );
    for (size_t j=0;j<valid-1;j++)
      for (size_t i=j+1;i<valid;i++)
	if ( area( bptr[j] ) > area( bptr[i] ) )
	  std::swap( bptr[j], bptr[i] );
  }

  void BVH4iBuilder::convertToSOALayout(const size_t threadID, const size_t numThreads)
  {
    const size_t startID = (threadID+0)*numNodes/numThreads;
    const size_t endID   = (threadID+1)*numNodes/numThreads;

    BVHNode  * __restrict__  bptr = ( BVHNode*)node + startID*4;

    BVH4i::Node * __restrict__  qptr = (BVH4i::Node*)node + startID;

    for (unsigned int n=startID;n<endID;n++,qptr++,bptr+=4)
      {
	prefetch<PFHINT_L1EX>(bptr+4);
	prefetch<PFHINT_L2EX>(bptr+4*4);
	convertToBVH4Layout(bptr);
	evictL1(bptr);
      }
  }

  void BVH4iBuilder::convertQBVHLayout(const size_t threadIndex, const size_t threadCount)
  {
    LockStepTaskScheduler::dispatchTask( task_convertToSOALayout, this, threadIndex, threadCount );    
  }

  __forceinline void computeAccelerationData(const unsigned int &geomID,
					     const unsigned int &primID,     
					     const Scene *__restrict__ const scene,
					     Triangle1 * __restrict__ const acc)
  {
    const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
    const TriangleMeshScene::TriangleMesh::Triangle & tri = mesh->triangle(primID);

    const mic_i pID(primID);
    const mic_i gID(geomID);

    const float *__restrict__ const vptr0 = (float*)&mesh->vertex(tri.v[0]);
    const float *__restrict__ const vptr1 = (float*)&mesh->vertex(tri.v[1]);
    const float *__restrict__ const vptr2 = (float*)&mesh->vertex(tri.v[2]);

    prefetch<PFHINT_L1>(vptr1);
    prefetch<PFHINT_L1>(vptr2);

    const mic_f v0 = broadcast4to16f(vptr0); //WARNING: zero last component
    const mic_f v1 = broadcast4to16f(vptr1);
    const mic_f v2 = broadcast4to16f(vptr2);

    const mic_f tri_accel = initTriangle1(v0,v1,v2,gID,pID,mic_i(mesh->mask));
    store16f_ngo(acc,tri_accel);
  }

  void BVH4iBuilder::createAccel(const size_t threadIndex, const size_t threadCount)
  {
    LockStepTaskScheduler::dispatchTask( task_createTriangle1Accel, this, threadIndex, threadCount );   
  }

  void BVH4iBuilder::createTriangle1Accel(const size_t threadID, const size_t numThreads)
  {
    const size_t startID = (threadID+0)*numPrimitives/numThreads;
    const size_t endID   = (threadID+1)*numPrimitives/numThreads;

    Triangle1    * __restrict__  acc  = accel + startID;
    const PrimRef* __restrict__  bptr = prims + startID;

    for (size_t j=startID; j<endID; j++, bptr++, acc++)
      {
	prefetch<PFHINT_NT>(bptr + L1_PREFETCH_ITEMS);
	prefetch<PFHINT_L2>(bptr + L2_PREFETCH_ITEMS);
	assert(bptr->geomID() < scene->size() );
	assert(bptr->primID() < scene->get( bptr->geomID() )->numPrimitives );

	computeAccelerationData(bptr->geomID(),bptr->primID(),scene,acc);
      }
  }


  void BVH4iBuilder::fillLocalWorkQueues(const size_t threadID, const size_t numThreads)
  {
    __aligned(64) BuildRecord br;
    const size_t numCores = (numThreads+3)/4;
    const size_t coreID   = threadID/4;

    if (threadID % 4 == 0)
      {
	unsigned int ID = coreID; // LockStepTaskScheduler::taskCounter.inc();
	
	while (true) 
	  {
	    /* get build record from global queue */
	    if (ID >= global_workStack.size()) break;
	    br = global_workStack.get(ID);
	    bool success = local_workStack[coreID].push(br);	  
	    if (!success) FATAL("can't fill local work queues");
	    ID += numCores;
	  }    

      }
  }

  void BVH4iBuilder::buildSubTrees(const size_t threadID, const size_t numThreads)
  {
    NodeAllocator alloc(atomicID,numAllocatedNodes);
    __aligned(64) BuildRecord br;
    const size_t numCores = (numThreads+3)/4;
    const size_t globalCoreID   = threadID/4;

#if defined(ENABLE_FILL_PER_CORE_WORK_QUEUES)

    if (numThreads > 1)
      {
	const size_t globalThreadID = threadID;
	const size_t localThreadID  = threadID % 4;
    
	if (localThreadID != 0)
	  {
	    localTaskScheduler[globalCoreID].dispatchTaskMainLoop(localThreadID,globalThreadID);
	  }
	else
	  {
	    local_workStack[globalCoreID].mutex.inc();
	    while (local_workStack[globalCoreID].size() < 8 && 
		   local_workStack[globalCoreID].size()+BVH4i::N <= SIZE_LOCAL_WORK_STACK) 
	      {
		BuildRecord br;
		if (!local_workStack[globalCoreID].pop_largest(br)) break;

		recurseSAH(br,alloc,FILL_LOCAL_QUEUES,globalThreadID,4);
	      }

	    localTaskScheduler[globalCoreID].releaseThreads(localThreadID,globalThreadID);	
	    local_workStack[globalCoreID].mutex.dec();
	  }
      }
#endif

    while(true)
      {
      /* process local work queue */
	while (1)
	  {
	    if (!local_workStack[globalCoreID].pop_largest(br)) 
	      {
		if (local_workStack[globalCoreID].mutex.val() > 0)
		  {
		    __pause(1024);
		    continue;
		  }
		else
		  break;
	      }
	    local_workStack[globalCoreID].mutex.inc();

	    recurseSAH(br,alloc,RECURSE,threadID,numThreads);
	    local_workStack[globalCoreID].mutex.dec();
	  }

	/* try task stealing */
        bool success = false;
#if defined(ENABLE_TASK_STEALING)
        for (size_t i=0; i<numThreads; i++)
        {
	  unsigned int next_threadID = (threadID+i);
	  if (next_threadID >= numThreads) next_threadID -= numThreads;
	  const unsigned int next_globalCoreID   = next_threadID/4;

	  assert(next_globalCoreID < numCores);
          if (local_workStack[next_globalCoreID].pop_smallest(br)) { 
            success = true;
            break;
          }
        }
#endif
        if (!success) break; 

	local_workStack[globalCoreID].mutex.inc();
	recurseSAH(br,alloc,RECURSE,threadID,numThreads);
	local_workStack[globalCoreID].mutex.dec();

      }

  }

  void reduceBinsParallel(const size_t currentThreadID,
			  const size_t childThreadID,
			  void *ptr)
  {
    Bin16 *__restrict__ bin16 = (Bin16*)ptr;
    bin16[childThreadID].prefetchL2();
    bin16[currentThreadID].merge(bin16[childThreadID]);
  }

  void BVH4iBuilder::parallelBinningGlobal(const size_t threadID, const size_t numThreads)
  {
    BuildRecord &current = global_sharedData.rec;

    const unsigned int items = current.items();
    const unsigned int startID = current.begin + ((threadID+0)*items/numThreads);
    const unsigned int endID   = current.begin + ((threadID+1)*items/numThreads);

    const mic_f centroidMin = broadcast4to16f(&current.bounds.centroid2.lower);
    const mic_f centroidMax = broadcast4to16f(&current.bounds.centroid2.upper);

    const mic_f centroidBoundsMin_2 = centroidMin;
    const mic_f centroidDiagonal_2  = centroidMax-centroidMin;
    const mic_f scale = select(centroidDiagonal_2 != 0.0f,rcp(centroidDiagonal_2) * mic_f(16.0f * 0.99f),mic_f::zero());

    PrimRef  *__restrict__ const tmp_prims = (PrimRef*)accel;

    fastbin_copy(prims,tmp_prims,startID,endID,centroidBoundsMin_2,scale,global_bin16[threadID]);    

    LockStepTaskScheduler::syncThreadsWithReduction( threadID, numThreads, reduceBinsParallel, global_bin16 );
    
    if (threadID == 0)
      {
	const float voxelArea = area(current.bounds.geometry);

	global_sharedData.split.cost = items * voxelArea;
	
	const Bin16 &bin16 = global_bin16[0];

	for (size_t dim=0;dim<3;dim++)
	  {
	    if (unlikely(centroidDiagonal_2[dim] == 0.0f)) continue;

	    const mic_f rArea = prefix_area_rl(bin16.min_x[dim],bin16.min_y[dim],bin16.min_z[dim],
					       bin16.max_x[dim],bin16.max_y[dim],bin16.max_z[dim]);
	    const mic_f lArea = prefix_area_lr(bin16.min_x[dim],bin16.min_y[dim],bin16.min_z[dim],
					       bin16.max_x[dim],bin16.max_y[dim],bin16.max_z[dim]);
	    const mic_i lnum  = prefix_count(bin16.count[dim]);

	    const mic_i rnum    = mic_i(items) - lnum;
	    const mic_i lblocks = (lnum + mic_i(3)) >> 2;
	    const mic_i rblocks = (rnum + mic_i(3)) >> 2;
	    const mic_m m_lnum  = lnum == 0;
	    const mic_m m_rnum  = rnum == 0;
	    const mic_f cost    = select(m_lnum|m_rnum,mic_f::inf(),lArea * mic_f(lblocks) + rArea * mic_f(rblocks) + voxelArea );

	    if (lt(cost,mic_f(global_sharedData.split.cost)))
	      {

		const mic_f min_cost    = vreduce_min(cost); 
		const mic_m m_pos       = min_cost == cost;
		const unsigned long pos = bitscan64(m_pos);	    
		
		assert(pos < 15);
		if (pos < 15)
		  {
		    global_sharedData.split.cost    = cost[pos];
		    global_sharedData.split.pos     = pos+1;
		    global_sharedData.split.dim     = dim;	    
		    global_sharedData.split.numLeft = lnum[pos];
		  }
	      }
	  }
      }
  }

  void BVH4iBuilder::parallelPartitioning(BuildRecord& current,
					  PrimRef * __restrict__ l_source,
					  PrimRef * __restrict__ r_source,
					  PrimRef * __restrict__ l_dest,
					  PrimRef * __restrict__ r_dest,
					  const Split &split,
					  Centroid_Scene_AABB &local_left,
					  Centroid_Scene_AABB &local_right)
  {
    const mic_f centroidMin = broadcast4to16f(&current.bounds.centroid2.lower);
    const mic_f centroidMax = broadcast4to16f(&current.bounds.centroid2.upper);

    const mic_f centroidBoundsMin_2 = centroidMin;
    const mic_f centroidDiagonal_2  = centroidMax-centroidMin;
    const mic_f scale = select(centroidDiagonal_2 != 0.0f,rcp(centroidDiagonal_2) * mic_f(16.0f * 0.99f),mic_f::zero());
 
    const unsigned int bestSplitDim = split.dim;
    const unsigned int bestSplit    = split.pos;

    const mic_f c = mic_f(centroidBoundsMin_2[bestSplitDim]);
    const mic_f s = mic_f(scale[bestSplitDim]);

    mic_f leftSceneBoundsMin((float)pos_inf);
    mic_f leftSceneBoundsMax((float)neg_inf);
    mic_f leftCentroidBoundsMin((float)pos_inf);
    mic_f leftCentroidBoundsMax((float)neg_inf);

    mic_f rightSceneBoundsMin((float)pos_inf);
    mic_f rightSceneBoundsMax((float)neg_inf);
    mic_f rightCentroidBoundsMin((float)pos_inf);
    mic_f rightCentroidBoundsMax((float)neg_inf);

    /* local queues for NGO stores */
    __aligned(64) PrimRef local_left_queue[2];
    __aligned(64) PrimRef local_right_queue[2];
    size_t num_local_left  = 0;
    size_t num_local_right = 0;

    const mic_m dim_mask = mic_m::shift1[bestSplitDim];

    for (;l_source<r_source;)
      {
	evictL1(l_source-2);	
	    
	prefetch<PFHINT_NT>(l_source+2);
	prefetch<PFHINT_L2>(l_source + L2_PREFETCH_ITEMS + 4);

	const mic_f b_min = broadcast4to16f(&l_source->lower);
	const mic_f b_max = broadcast4to16f(&l_source->upper);
	const mic_f b_centroid2 = b_min + b_max;

	if (likely(lt_split(b_min,b_max,dim_mask,c,s,mic_f(bestSplit)))) 
	  {
	    store4f(&local_left_queue[num_local_left].lower,b_min);
	    store4f(&local_left_queue[num_local_left].upper,b_max);
	    num_local_left++;
	    l_source++;

	    leftSceneBoundsMin = min(leftSceneBoundsMin,b_min);
	    leftSceneBoundsMax = max(leftSceneBoundsMax,b_max);

	    leftCentroidBoundsMin = min(leftCentroidBoundsMin,b_centroid2);
	    leftCentroidBoundsMax = max(leftCentroidBoundsMax,b_centroid2);
	    
	    if (unlikely(((size_t)l_dest % 64) != 0) && num_local_left == 1) 
	      {
		*l_dest++ = local_left_queue[0];
		num_local_left--;
	      }
	    else
	      {
		if (num_local_left == 2)
		  {
		    const mic_f twoAABBs = load16f(local_left_queue);
		    num_local_left = 0;
		    store16f_ngo(l_dest,twoAABBs);
		    l_dest+=2;
		  }
	      }	

	  }
	else
	  {
	    store4f(&local_right_queue[num_local_right].lower,b_min);
	    store4f(&local_right_queue[num_local_right].upper,b_max);
	    num_local_right++;
	    l_source++;

	    rightSceneBoundsMin = min(rightSceneBoundsMin,b_min);
	    rightSceneBoundsMax = max(rightSceneBoundsMax,b_max);
	    rightCentroidBoundsMin = min(rightCentroidBoundsMin,b_centroid2);
	    rightCentroidBoundsMax = max(rightCentroidBoundsMax,b_centroid2);
	    
	    if (unlikely(((size_t)r_dest % 64) != 0) && num_local_right == 1)
	      {
		*r_dest++ = local_right_queue[0];
		num_local_right--;
	      }
	    else
	      {
		if (num_local_right == 2)
		  {
		    const mic_f twoAABBs = load16f(local_right_queue);
		    num_local_right = 0;
		    store16f_ngo(r_dest,twoAABBs);
		    r_dest+=2;
		  }
	      }	
	  }
      }

    /* flush local queues */

    if (num_local_left % 2 != 0)
      *l_dest = local_left_queue[0];

    if (num_local_right % 2 != 0)
      *r_dest = local_right_queue[0];
    
    store4f(&local_left.geometry.lower,leftSceneBoundsMin);
    store4f(&local_left.geometry.upper,leftSceneBoundsMax);
    store4f(&local_left.centroid2.lower,leftCentroidBoundsMin);
    store4f(&local_left.centroid2.upper,leftCentroidBoundsMax);

    store4f(&local_right.geometry.lower,rightSceneBoundsMin);
    store4f(&local_right.geometry.upper,rightSceneBoundsMax);
    store4f(&local_right.centroid2.lower,rightCentroidBoundsMin);
    store4f(&local_right.centroid2.upper,rightCentroidBoundsMax);


  }

  void BVH4iBuilder::parallelPartitioningGlobal(const size_t threadID, const size_t numThreads)
  {
    BuildRecord &current = global_sharedData.rec;

    const unsigned int items = current.items();
    const unsigned int startID = current.begin + ((threadID+0)*items/numThreads);
    const unsigned int endID   = current.begin + ((threadID+1)*items/numThreads);
   
 
    const unsigned int bestSplitDim     = global_sharedData.split.dim;
    const unsigned int bestSplit        = global_sharedData.split.pos;
    const unsigned int bestSplitNumLeft = global_sharedData.split.numLeft;


    const mic_i lnum    = prefix_sum(global_bin16[threadID].thread_count[bestSplitDim]);
    const unsigned int local_numLeft = lnum[bestSplit-1];
    const unsigned int local_numRight = (endID-startID) - lnum[bestSplit-1];
 
    const unsigned int thread_start_left  = global_sharedData.lCounter.add(local_numLeft);
    const unsigned int thread_start_right = global_sharedData.rCounter.add(local_numRight);

    PrimRef  *__restrict__ const tmp_prims = (PrimRef*)accel;

    PrimRef * __restrict__ l_source = tmp_prims + startID;
    PrimRef * __restrict__ r_source = tmp_prims + endID;

    PrimRef * __restrict__ l_dest     = prims + current.begin + thread_start_left;
    PrimRef * __restrict__ r_dest     = prims + current.begin + thread_start_right + bestSplitNumLeft;

    __aligned(64) Centroid_Scene_AABB local_left;
    __aligned(64) Centroid_Scene_AABB local_right; // just one local

    parallelPartitioning(current,l_source,r_source,l_dest,r_dest,global_sharedData.split,local_left,local_right);

    global_sharedData.left.extend_atomic(local_left); 
    global_sharedData.right.extend_atomic(local_right);  
  }


  // =======================================================================================================
  // =======================================================================================================
  // =======================================================================================================

  bool BVH4iBuilder::splitSequential(BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild)
  {
#if defined(DEBUG)
    checkBuildRecord(current);
#endif

    /* mark as leaf if leaf threshold reached */
    if (current.items() <= BVH4i::N) {
      current.createLeaf();
      return false;
    }
    
    const mic_f centroidMin = broadcast4to16f(&current.bounds.centroid2.lower);
    const mic_f centroidMax = broadcast4to16f(&current.bounds.centroid2.upper);

    const mic_f centroidBoundsMin_2 = centroidMin;
    const mic_f centroidDiagonal_2  = centroidMax-centroidMin;
    const mic_f scale = select(centroidDiagonal_2 != 0.0f,rcp(centroidDiagonal_2) * mic_f(16.0f * 0.99f),mic_f::zero());

    mic_f leftArea[3];
    mic_f rightArea[3];
    mic_i leftNum[3];

    fastbin(prims,current.begin,current.end,centroidBoundsMin_2,scale,leftArea,rightArea,leftNum);

    const unsigned int items = current.items();
    const float voxelArea = area(current.bounds.geometry);
    Split split;
    split.cost = items * voxelArea;

    for (size_t dim = 0;dim < 3;dim++) 
      {
	if (unlikely(centroidDiagonal_2[dim] == 0.0f)) continue;

	const mic_f rArea   = rightArea[dim]; // bin16.prefix_area_rl(dim);
	const mic_f lArea   = leftArea[dim];  // bin16.prefix_area_lr(dim);      
	const mic_i lnum    = leftNum[dim];   // bin16.prefix_count(dim);

	const mic_i rnum    = mic_i(items) - lnum;
	const mic_i lblocks = (lnum + mic_i(3)) >> 2;
	const mic_i rblocks = (rnum + mic_i(3)) >> 2;
	const mic_m m_lnum  = lnum == 0;
	const mic_m m_rnum  = rnum == 0;
	const mic_f cost    = select(m_lnum|m_rnum,mic_f::inf(),lArea * mic_f(lblocks) + rArea * mic_f(rblocks) + voxelArea );

	if (lt(cost,mic_f(split.cost)))
	  {

	    const mic_f min_cost    = vreduce_min(cost); 
	    const mic_m m_pos       = min_cost == cost;
	    const unsigned long pos = bitscan64(m_pos);	    

	    assert(pos < 15);

	    if (pos < 15)
	      {
		split.cost    = cost[pos];
		split.pos     = pos+1;
		split.dim     = dim;	    
		split.numLeft = lnum[pos];
	      }
	  }
      };

    if (unlikely(split.pos == -1)) 
      split_fallback(prims,current,leftChild,rightChild);
   // /* partitioning of items */
    else 
      {
	leftChild.bounds.reset();
	rightChild.bounds.reset();

	const unsigned int mid = partitionPrimRefs<L2_PREFETCH_ITEMS>(prims ,current.begin, current.end-1, split.pos, split.dim, centroidBoundsMin_2, scale, leftChild.bounds, rightChild.bounds);

	assert(area(leftChild.bounds.geometry) >= 0.0f);

#if defined(DEBUG)
	if (current.begin + mid != current.begin + split.numLeft)
	  {
	    mtx.lock();	    
	    DBG_PRINT(current);
	    DBG_PRINT(mid);
	    DBG_PRINT(split);
	    DBG_PRINT(leftNum[0]);
	    DBG_PRINT(leftNum[1]);
	    DBG_PRINT(leftNum[2]);

	    checkBuildRecord(current);
	    
	    mtx.unlock();
	  }
#endif

	assert(current.begin + mid == current.begin + split.numLeft);

	if (unlikely(current.begin + mid == current.begin || current.begin + mid == current.end)) 
	  {
	    std::cout << "WARNING: mid == current.begin || mid == current.end " << std::endl;
	    DBG_PRINT(split);
	    DBG_PRINT(current);
	    DBG_PRINT(mid);
	    split_fallback(prims,current,leftChild,rightChild);	    
	  }
	else
	  {
	    const unsigned int current_mid = current.begin + split.numLeft;
	    leftChild.init(current.begin,current_mid);
	    rightChild.init(current_mid,current.end);
	  }

      }


#if defined(DEBUG)
    checkBuildRecord(leftChild);
    checkBuildRecord(rightChild);
#endif

    if (leftChild.items()  <= BVH4i::N) leftChild.createLeaf();
    if (rightChild.items() <= BVH4i::N) rightChild.createLeaf();	
    return true;
  }

  bool BVH4iBuilder::splitParallelGlobal( BuildRecord &current,
					  BuildRecord &leftChild,
					  BuildRecord &rightChild,
					  const size_t threadID,
					  const size_t numThreads)
  {
    const unsigned int items = current.end - current.begin;
    assert(items >= BUILD_RECORD_PARALLEL_SPLIT_THRESHOLD);

#if defined(DEBUG)
    checkBuildRecord(current);
#endif
  
    /* mark as leaf if leaf threshold reached */
    if (items <= BVH4i::N) {
      current.createLeaf();
      return false;
    }

     global_sharedData.rec = current;
     global_sharedData.split.reset();
     global_sharedData.left.reset();
     global_sharedData.right.reset();
     
     LockStepTaskScheduler::dispatchTask( task_parallelBinningGlobal, this, threadID, numThreads );

     if (unlikely(global_sharedData.split.pos == -1)) 
       split_fallback(prims,current,leftChild,rightChild);
     else
       {
	 global_sharedData.left.reset();
	 global_sharedData.right.reset();

	 global_sharedData.lCounter.reset(0);
	 global_sharedData.rCounter.reset(0); 

	 LockStepTaskScheduler::dispatchTask( task_parallelPartitioningGlobal, this, threadID, numThreads );

	 const unsigned int mid = current.begin + global_sharedData.split.numLeft;

	if (unlikely(current.begin == mid || mid == current.end)) 
	  {
	    std::cout << "WARNING: mid == current.begin || mid == current.end " << std::endl;
	    DBG_PRINT(global_sharedData.split);
	    DBG_PRINT(current);
	    DBG_PRINT(mid);
	    split_fallback(prims,current,leftChild,rightChild);	    
	  }
	else
	  {
	    leftChild.init(global_sharedData.left,current.begin,mid);
	    rightChild.init(global_sharedData.right,mid,current.end);
	  }	 
       }

#if defined(DEBUG)
     checkBuildRecord(leftChild);
     checkBuildRecord(rightChild);
#endif
     
     if (leftChild.items()  <= BVH4i::N) leftChild.createLeaf();
     if (rightChild.items() <= BVH4i::N) rightChild.createLeaf();
     return true;
  }



  bool BVH4iBuilder::splitParallelLocal(BuildRecord &current,
					BuildRecord &leftChild,
					BuildRecord &rightChild,
					const size_t threadID)
  {
    const unsigned int items    = current.end - current.begin;
    const size_t globalCoreID   = threadID / 4;
    const size_t localThreadID  = threadID % 4;
    const size_t globalThreadID = threadID;
    
    assert(items >= THRESHOLD_FOR_SUBTREE_RECURSION);

#if defined(DEBUG)
    checkBuildRecord(current);
#endif
  
    /* mark as leaf if leaf threshold reached */
    if (items <= BVH4i::N) {
      current.createLeaf();
      return false;
    }

    SharedBinningPartitionData &sd = local_sharedData[globalCoreID]; 

    sd.rec = current;
    sd.split.reset();
    sd.left.reset();
    sd.right.reset();

    localTaskScheduler[globalCoreID].dispatchTask( task_parallelBinningLocal, this, localThreadID, globalThreadID );

    if (unlikely(sd.split.pos == -1)) 
      split_fallback(prims,current,leftChild,rightChild);
    else
      {

	 sd.left.reset();
	 sd.right.reset();

	 sd.lCounter.reset(0);
	 sd.rCounter.reset(0); 

	 localTaskScheduler[globalCoreID].dispatchTask( task_parallelPartitioningLocal, this, localThreadID, globalThreadID );

	 const unsigned int mid = current.begin + sd.split.numLeft;

	 if (unlikely(mid == current.begin || mid == current.end)) 
	   {
	     std::cout << "WARNING: mid == current.begin || mid == current.end " << std::endl;
	     DBG_PRINT(sd.split);
	     DBG_PRINT(current);
	     DBG_PRINT(mid);
	     split_fallback(prims,current,leftChild,rightChild);	    
	   }
	 else
	   {
	     leftChild.init(sd.left,current.begin,mid);
	     rightChild.init(sd.right,mid,current.end);
	   }
	 
       }

#if defined(DEBUG)
    checkBuildRecord(leftChild);
    checkBuildRecord(rightChild);
#endif
     
     if (leftChild.items()  <= BVH4i::N) leftChild.createLeaf();
     if (rightChild.items() <= BVH4i::N) rightChild.createLeaf();
     return true;
  }


  __forceinline bool BVH4iBuilder::split(BuildRecord& current, BuildRecord& left, BuildRecord& right, const size_t mode, const size_t threadID, const size_t numThreads)
  {
    if (unlikely(mode == BUILD_TOP_LEVEL))
      {
	if (current.items() >= BUILD_RECORD_PARALLEL_SPLIT_THRESHOLD)
	  return splitParallelGlobal(current,left,right,threadID,numThreads);
	else
	  {
	    DBG(std::cout << "WARNING in top-level build: too few items for parallel split " << current.items() << std::endl << std::flush);
	    return splitSequential(current,left,right);
	  }
      }
    else if (unlikely(mode == FILL_LOCAL_QUEUES))
      {
	if (current.items() >= THRESHOLD_FOR_SUBTREE_RECURSION)
	  return splitParallelLocal(current,left,right,threadID);
	else
	  {
	    DBG(std::cout << "WARNING in fill_local_queues build: too few items for parallel split " << current.items() << std::endl << std::flush);
	    return splitSequential(current,left,right);
	  }
	
      }
    else
      return splitSequential(current,left,right);
  }

  // =======================================================================================================
  // =======================================================================================================
  // =======================================================================================================
  
  __forceinline void BVH4iBuilder::createLeaf(BuildRecord& current, NodeAllocator& alloc,const size_t threadIndex, const size_t threadCount)
  {
#if defined(DEBUG)
    if (current.depth > BVH4i::maxBuildDepthLeaf) 
      throw std::runtime_error("ERROR: depth limit reached");
#endif
    
    /* create leaf */
    if (current.items() <= BVH4i::N) {
      node[current.parentID].createLeaf(current.begin,current.items());
#if defined(DEBUG)
      checkLeafNode(node[current.parentID]);      
#endif
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
    const size_t currentIndex = alloc.get(numNodesToAllocate);

    node[current.parentID].createNode(currentIndex,numChildren);
    
    /* recurse into each child */
    for (size_t i=0; i<numChildren; i++) 
    {
      node[currentIndex+i].lower = children[i].bounds.geometry.lower;
      node[currentIndex+i].upper = children[i].bounds.geometry.upper;
      children[i].parentID = currentIndex+i;
      children[i].depth = current.depth+1;
      createLeaf(children[i],alloc,threadIndex,threadCount);
    }
  }  

  __forceinline void BVH4iBuilder::recurse(BuildRecord& current, NodeAllocator& alloc,const size_t mode, const size_t threadID, const size_t numThreads)
  {
    if (mode == BUILD_TOP_LEVEL) {
      global_workStack.push_nolock(current);
    }
    else if (current.items() > THRESHOLD_FOR_SUBTREE_RECURSION || mode == FILL_LOCAL_QUEUES) {
      const size_t coreID = threadID/4;
      if (!local_workStack[coreID].push(current))
        recurseSAH(current,alloc,RECURSE,threadID,numThreads);
    }
    else
      recurseSAH(current,alloc,RECURSE,threadID,numThreads);
  }
  
  void BVH4iBuilder::recurseSAH(BuildRecord& current, NodeAllocator& alloc,const size_t mode, const size_t threadID, const size_t numThreads)
  {
    __aligned(64) BuildRecord children[BVH4i::N];

    /* create leaf node */
    if (current.depth >= BVH4i::maxBuildDepth || current.isLeaf()) {
      createLeaf(current,alloc,threadID,numThreads);
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
      createLeaf(current,alloc,threadID,numThreads);
      return;
    }

    /* allocate next four nodes */
    const size_t currentIndex = alloc.get(numNodesToAllocate);
    node[current.parentID].createNode(currentIndex,numChildren);

    /* init used/unused nodes */
    const mic_f init_node = load16f((float*)BVH4i::initQBVHNode);
    store16f_ngo((float*)&node[currentIndex+0],init_node);
    store16f_ngo((float*)&node[currentIndex+2],init_node);

    /* recurse into each child */
    for (unsigned int i=0; i<numChildren; i++) 
    {
      node[currentIndex+i].lower = (Vec3fa) children[i].bounds.geometry.lower;
      node[currentIndex+i].upper = (Vec3fa) children[i].bounds.geometry.upper;
      children[i].parentID = currentIndex+i;
      recurse(children[i],alloc,mode,threadID,numThreads);
    }

  }

  void BVH4iBuilder::checkLeafNode(const BVHNode &entry)
  {
    if (!entry.isLeaf())
      FATAL("no leaf");

    unsigned int accel_entries = entry.items();
    unsigned int accel_offset  = entry.itemListOfs();

    BBox3f leaf_prim_bounds = empty;
    for (size_t i=0;i<accel_entries;i++)
      {
	leaf_prim_bounds.extend( prims[ accel_offset + i ].lower );
	leaf_prim_bounds.extend( prims[ accel_offset + i ].upper );
      }

    BBox3f leaf_tri_bounds = empty;
    for (size_t i=0;i<accel_entries;i++)
      {
	const unsigned int geomID = prims[ accel_offset + i ].geomID();
	const unsigned int primID = prims[ accel_offset + i ].primID();

	const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
	const TriangleMeshScene::TriangleMesh::Triangle & tri = mesh->triangle(primID);

	leaf_tri_bounds.extend( mesh->vertex(tri.v[0]) );
	leaf_tri_bounds.extend( mesh->vertex(tri.v[1]) );
	leaf_tri_bounds.extend( mesh->vertex(tri.v[2]) );	
      }

    if (!(subset(leaf_prim_bounds,leaf_tri_bounds) && subset(leaf_tri_bounds,leaf_prim_bounds))) 
      {
	DBG_PRINT(leaf_prim_bounds);
	DBG_PRINT(leaf_tri_bounds);
	FATAL("check build record");
      }

  }


  void BVH4iBuilder::checkBuildRecord(const BuildRecord &current)
  {
    BBox3f check_box;
    BBox3f box;
    check_box = empty;
    box = *(BBox3f*)&current.bounds.geometry;

    BBox3f *aabb = (BBox3f*)prims;

    for (unsigned int i=current.begin;i<current.end;i++) 
      {
	check_box.extend(aabb[i]);
	if (!subset(aabb[i],box))
	  {
	    DBG_PRINT(current);
	    DBG_PRINT(i);
	    DBG_PRINT(prims[i]);
	    FATAL("check build record => subset");
	  }
      }

    //if (enablePreSplits) return;
    if (!(subset(check_box,box) && subset(box,check_box))) 
      {
	DBG_PRINT(current);
	DBG_PRINT(check_box);
	DBG_PRINT(box);
	FATAL("check build record");
      }
  }

  void BVH4iBuilder::parallelBinningLocal(const size_t localThreadID,const size_t globalThreadID)
  {
    const size_t globalCoreID = globalThreadID/4;
    BuildRecord &current = local_sharedData[globalCoreID].rec;

    const unsigned int items   = current.items();
    const unsigned int startID = current.begin + ((localThreadID+0)*items/4);
    const unsigned int endID   = current.begin + ((localThreadID+1)*items/4);
    
    const mic_f centroidMin = broadcast4to16f(&current.bounds.centroid2.lower);
    const mic_f centroidMax = broadcast4to16f(&current.bounds.centroid2.upper);

    const mic_f centroidBoundsMin_2 = centroidMin;
    const mic_f centroidDiagonal_2  = centroidMax-centroidMin;
    const mic_f scale = select(centroidDiagonal_2 != 0.0f,rcp(centroidDiagonal_2) * mic_f(16.0f * 0.99f),mic_f::zero());

    PrimRef  *__restrict__ const tmp_prims = (PrimRef*)accel;

    fastbin_copy(prims,tmp_prims,startID,endID,centroidBoundsMin_2,scale,global_bin16[globalThreadID]);    

    localTaskScheduler[globalCoreID].syncThreads(localThreadID);

    if (localThreadID == 0)
      {
	Bin16 &bin16 = global_bin16[globalThreadID];

	for (size_t i=1;i<4;i++)
	  bin16.merge(global_bin16[globalThreadID+i]);

	const float voxelArea = area(current.bounds.geometry);

	local_sharedData[globalCoreID].split.cost = items * voxelArea;	

	for (size_t dim=0;dim<3;dim++)
	  {
	    if (unlikely(centroidDiagonal_2[dim] == 0.0f)) continue;

	    const mic_f rArea = prefix_area_rl(bin16.min_x[dim],bin16.min_y[dim],bin16.min_z[dim],
					       bin16.max_x[dim],bin16.max_y[dim],bin16.max_z[dim]);
	    const mic_f lArea = prefix_area_lr(bin16.min_x[dim],bin16.min_y[dim],bin16.min_z[dim],
					       bin16.max_x[dim],bin16.max_y[dim],bin16.max_z[dim]);
	    const mic_i lnum  = prefix_count(bin16.count[dim]);

	    const mic_i rnum    = mic_i(items) - lnum;
	    const mic_i lblocks = (lnum + mic_i(3)) >> 2;
	    const mic_i rblocks = (rnum + mic_i(3)) >> 2;
	    const mic_m m_lnum  = lnum == 0;
	    const mic_m m_rnum  = rnum == 0;
	    const mic_f cost    = select(m_lnum|m_rnum,mic_f::inf(),lArea * mic_f(lblocks) + rArea * mic_f(rblocks) + voxelArea );

	    if (lt(cost,mic_f(local_sharedData[globalCoreID].split.cost)))
	      {

		const mic_f min_cost    = vreduce_min(cost); 
		const mic_m m_pos       = min_cost == cost;
		const unsigned long pos = bitscan64(m_pos);	    

		assert(pos < 15);
		if (pos < 15)
		  {
		    local_sharedData[globalCoreID].split.cost    = cost[pos];
		    local_sharedData[globalCoreID].split.pos     = pos+1;
		    local_sharedData[globalCoreID].split.dim     = dim;	    
		    local_sharedData[globalCoreID].split.numLeft = lnum[pos];
		  }
	      }
	  }
      }

  }


  void BVH4iBuilder::parallelPartitioningLocal(const size_t localThreadID,const size_t globalThreadID)
  {
    const size_t threads = 4;
    const size_t globalCoreID = globalThreadID/threads;

    SharedBinningPartitionData &sd = local_sharedData[globalCoreID];    
    BuildRecord &current = sd.rec;
    Bin16 &bin16 = global_bin16[globalThreadID];

    // ----------------------------------------------

    const unsigned int items = current.items();
    const unsigned int startID = current.begin + ((localThreadID+0)*items/threads);
    const unsigned int endID   = current.begin + ((localThreadID+1)*items/threads);
   
    const unsigned int bestSplitDim     = sd.split.dim;
    const unsigned int bestSplit        = sd.split.pos;
    const unsigned int bestSplitNumLeft = sd.split.numLeft;

    const mic_i lnum    = prefix_sum(bin16.thread_count[bestSplitDim]);
    const unsigned int local_numLeft = lnum[bestSplit-1];
    const unsigned int local_numRight = (endID-startID) - lnum[bestSplit-1];
 
    const unsigned int thread_start_left  = sd.lCounter.add(local_numLeft);
    const unsigned int thread_start_right = sd.rCounter.add(local_numRight);

    PrimRef  *__restrict__ const tmp_prims = (PrimRef*)accel;

    PrimRef * __restrict__ l_source = tmp_prims + startID;
    PrimRef * __restrict__ r_source = tmp_prims + endID;

    PrimRef * __restrict__ l_dest     = prims + current.begin + thread_start_left;
    PrimRef * __restrict__ r_dest     = prims + current.begin + thread_start_right + bestSplitNumLeft;

    __aligned(64) Centroid_Scene_AABB local_left;
    __aligned(64) Centroid_Scene_AABB local_right; 

    parallelPartitioning(current,l_source,r_source,l_dest,r_dest,sd.split,local_left,local_right);

    sd.left.extend_atomic(local_left); 
    sd.right.extend_atomic(local_right);  
    
  }

  void BVH4iBuilder::computePrimRefs(const size_t threadIndex, const size_t threadCount)
  {
    LockStepTaskScheduler::dispatchTask( task_computePrimRefsTriangles, this, threadIndex, threadCount );
  }


  // =======================================================================================================
  // =======================================================================================================
  // =======================================================================================================

  bool split_fallback(PrimRef * __restrict__ const primref, BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild)
  {
    const unsigned int center = (current.begin + current.end)/2;
    
    Centroid_Scene_AABB left; left.reset();
    for (size_t i=current.begin; i<center; i++)
      left.extend(primref[i].bounds());
    leftChild.init(left,current.begin,center);
    
    Centroid_Scene_AABB right; right.reset();
    for (size_t i=center; i<current.end; i++)
      right.extend(primref[i].bounds());	
    rightChild.init(right,center,current.end);
    
    return true;
  }


  void BVH4iBuilder::build_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event) 
  {
    TIMER(double msec = 0.0);

    /* initialize thread-local work stacks */
    if (threadIndex % 4 == 0)
      local_workStack[threadIndex].reset();

    /* all worker threads enter tasking system */
    if (threadIndex != 0) {
      LockStepTaskScheduler::dispatchTaskMainLoop(threadIndex,threadCount); 
      return;
    }

    /* start measurement */
    double t0 = 0.0f;
#if !defined(PROFILE)
    if (g_verbose >= 2) 
#endif
      t0 = getSeconds();

    TIMER(msec = getSeconds());
    
    /* calculate list of primrefs */
    global_bounds.reset();

    computePrimRefs(threadIndex,threadCount);

    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "task_computePrimRefs " << 1000. * msec << " ms" << std::endl << std::flush);
    TIMER(msec = getSeconds());

    /* allocate and initialize root node */
    atomicID.reset(numNodesToAllocate);
    node[0].lower = global_bounds.geometry.lower;
    node[0].upper = global_bounds.geometry.upper;
    
    /* create initial build record */
    BuildRecord br;
    br.init(global_bounds,0,numPrimitives);
    br.depth = 1;
    br.parentID = 0;
        
    /* push initial build record to global work stack */
    global_workStack.reset();
    global_workStack.push_nolock(br);    

    /* work in multithreaded toplevel mode until sufficient subtasks got generated */    
    NodeAllocator alloc(atomicID,numAllocatedNodes);
    const size_t coreCount = (threadCount+3)/4;
    while (global_workStack.size() < coreCount &&
	   global_workStack.size()+BVH4i::N <= SIZE_GLOBAL_WORK_STACK) 
    {
      BuildRecord br;
      if (!global_workStack.pop_nolock_largest(br)) break;
      recurseSAH(br,alloc,BUILD_TOP_LEVEL,threadIndex,threadCount);
    }

    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "build_top_level " << 1000. * msec << " ms" << std::endl << std::flush);

    /* fill per core work queues */    
    TIMER(msec = getSeconds());    
    LockStepTaskScheduler::dispatchTask(task_fillLocalWorkQueues, this, threadIndex, threadCount );
    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "task_fillLocalWorkQueues " << 1000. * msec << " ms" << std::endl << std::flush);

    /* now process all created subtasks on multiple threads */    
    TIMER(msec = getSeconds());    
    LockStepTaskScheduler::dispatchTask(task_buildSubTrees, this, threadIndex, threadCount );
    numNodes = atomicID >> 2;
    DBG(DBG_PRINT(atomicID));
    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "task_buildSubTrees " << 1000. * msec << " ms" << std::endl << std::flush);

    /* create triangle acceleration structure */
    TIMER(msec = getSeconds());        
    createAccel(threadIndex, threadCount );
    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "task_createAccel " << 1000. * msec << " ms" << std::endl << std::flush);

    
    /* convert to SOA node layout */
    TIMER(msec = getSeconds());     
    convertQBVHLayout(threadIndex, threadCount );
    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "task_convertToSOALayout " << 1000. * msec << " ms" << std::endl << std::flush);

    
    /* update BVH4 */
    bvh->root = bvh->qbvh[0].lower[0].child; 
    bvh->bounds = BBox3f(*(Vec3fa*)&bvh->qbvh->lower[0],*(Vec3fa*)&bvh->qbvh->upper[0]);
    
    /* release all threads again */
    LockStepTaskScheduler::releaseThreads(threadCount);

    /* stop measurement */
#if !defined(PROFILE)
    if (g_verbose >= 2) 
#endif
      dt = getSeconds()-t0;
  }
};
