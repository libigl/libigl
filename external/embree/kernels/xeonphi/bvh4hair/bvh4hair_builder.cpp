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

#include "bvh4hair/bvh4hair_builder.h"
#include "geometry/bezier1i.h"
#include "bvh4hair/bvh4hair_statistics.h"

namespace embree
{
#define DBG(x) 

#define L1_PREFETCH_ITEMS 4
#define L2_PREFETCH_ITEMS 16
#define SINGLE_THREADED_BUILD_THRESHOLD        512
#define THRESHOLD_FOR_SUBTREE_RECURSION         64
#define BUILD_RECORD_PARALLEL_SPLIT_THRESHOLD 1024

#define ENABLE_OBB_BVH4 1
#define ENABLE_AABB_NODES 1

#define BVH4HAIR_NODE_PREALLOC_FACTOR  0.7f
#define INTERSECTION_COST              1.0f

#define AABB_OBB_SWITCH_THRESHOLD      1.1f
#define THRESHOLD_SWITCH

#define TIMER(x)  

#if defined(DEBUG)
  extern AtomicMutex mtx;
#endif

  static double dt = 0.0f;

  // ==========================================================================================
  // ==========================================================================================
  // ==========================================================================================

  enum {
    BVH4HAIR_AABB_NODE = 1,
    BVH4HAIR_OBB_NODE  = 2
  };

  void BVH4HairBuilder::printBuilderName()
  {
    std::cout << "building BVH4Hair using binned SAH builder (MIC) ... " << std::endl;    
    
  }

  size_t BVH4HairBuilder::getNumPrimitives()
  {
    DBG(PING);

    /* count total number of virtual objects */
    size_t numCurves = 0;       
    for (size_t i=0;i<scene->size();i++)
      {
	if (unlikely(scene->get(i) == NULL)) continue;
	if (unlikely((scene->get(i)->type != BEZIER_CURVES))) continue;
	if (unlikely(!scene->get(i)->isEnabled())) continue;
        BezierCurves* geom = (BezierCurves*) scene->getBezierCurves(i);
	numCurves += geom->numCurves;
      }
    return numCurves;	
  }

  void BVH4HairBuilder::allocateData(const size_t threadCount, const size_t totalNumPrimitives)
  {
    DBG(PING);
    size_t numPrimitivesOld = numPrimitives;
    numPrimitives = totalNumPrimitives;

    if (numPrimitivesOld != numPrimitives)
      {
	const size_t numPrims = numPrimitives;
	const size_t minAllocNodes = (threadCount+1) * ALLOCATOR_NODE_BLOCK_SIZE;
	size_t numNodes = (size_t)((numPrims+1)/2 * BVH4HAIR_NODE_PREALLOC_FACTOR) + minAllocNodes;
	if (numPrimitives == 0) numNodes = 0;
	allocateMemoryPools(numPrims,numNodes);
      }
  }

  void BVH4HairBuilder::computePrimRefs(const size_t threadIndex, const size_t threadCount)
  {
    DBG(PING);
    scene->lockstep_scheduler.dispatchTask( task_computePrimRefsBezierCurves, this, threadIndex, threadCount );	
  }

  void BVH4HairBuilder::computePrimRefsBezierCurves(const size_t threadID, const size_t numThreads) 
  {
    DBG(PING);

    const size_t numTotalGroups = scene->size();

    /* count total number of virtual objects */
    const size_t numBezierCurves = numPrimitives;
    const size_t startID   = (threadID+0)*numBezierCurves/numThreads;
    const size_t endID     = (threadID+1)*numBezierCurves/numThreads; 
        
    Bezier1i *__restrict__ const bptr     = (Bezier1i*)this->prims;

    // === find first group containing startID ===
    unsigned int g=0, numSkipped = 0;
    for (; g<numTotalGroups; g++) {       
      if (unlikely(scene->get(g) == NULL)) continue;
      if (unlikely((scene->get(g)->type != BEZIER_CURVES))) continue;
      if (unlikely(!scene->get(g)->isEnabled())) continue;
      BezierCurves* geom = (BezierCurves*) scene->getBezierCurves(g);
      const size_t numPrims = geom->numCurves;
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
	if (unlikely((scene->get(g)->type != BEZIER_CURVES))) continue;
	if (unlikely(!scene->get(g)->isEnabled())) continue;

	BezierCurves* geom = (BezierCurves*) scene->getBezierCurves(g);

        size_t N = geom->numCurves;
        for (unsigned int i=offset; i<N && currentID < endID; i++, currentID++)	 
        { 			    
	  const mic2f b2 = geom->bounds_mic2f(i);
	  const mic_f bmin = b2.x;
	  const mic_f bmax = b2.y;
          
          bounds_scene_min = min(bounds_scene_min,bmin);
          bounds_scene_max = max(bounds_scene_max,bmax);
          const mic_f centroid2 = bmin+bmax;
          bounds_centroid_min = min(bounds_centroid_min,centroid2);
          bounds_centroid_max = max(bounds_centroid_max,centroid2);

	  bptr[currentID].p = geom->fristVertexPtr(i); // FIXME: this does not support strides!!
          bptr[currentID].geomID = g;
          bptr[currentID].primID = i;
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


  void BVH4HairBuilder::allocateMemoryPools(const size_t numPrims,
					    const size_t numNodes)
  {
    const size_t additional_size = 16 * CACHELINE_SIZE;
    const size_t sizeNodeInBytes    = sizeof(BVH4Hair::UnalignedNode);
    const size_t sizePrimRefInBytes = sizeof(Bezier1i);
    const size_t sizeAccelInBytes   = sizeof(Bezier1i);

    /* free previously allocated memory */

    if (prims)  {
      assert(size_prims > 0);
      os_free(prims,size_prims);
    }
    if (node  ) {
      assert(bvh4hair->size_node > 0);
      os_free(node ,bvh4hair->size_node);
    }
    if (accel ) {
      assert(bvh4hair->size_accel > 0);
      os_free(accel,bvh4hair->size_accel);
    }
      
    // === allocated memory for primrefs,nodes, and accel ===
    const size_t size_primrefs = numPrims * sizePrimRefInBytes + additional_size;
    const size_t size_node     = (double)(numNodes * sizeNodeInBytes + additional_size) * g_memory_preallocation_factor;
    const size_t size_accel    = numPrims * sizeAccelInBytes   + additional_size;

    numAllocated64BytesBlocks = size_node / sizeof(BVH4Hair::UnalignedNode); // FIXME: do memory handling in 64 byte blocks
      
#if DEBUG
    DBG_PRINT(numAllocated64BytesBlocks);
    DBG_PRINT(sizeNodeInBytes);
    DBG_PRINT(sizePrimRefInBytes);
    DBG_PRINT(sizeAccelInBytes);

    DBG_PRINT(size_primrefs);
    DBG_PRINT(size_node);
    DBG_PRINT(size_accel);
#endif

    prims = (Bezier1i                 *) os_malloc(size_primrefs); 
    node  = (BVH4Hair::UnalignedNode  *) os_malloc(size_node);
    accel = (Bezier1i                 *) os_malloc(size_accel);

    assert(prims  != 0);
    assert(node   != 0);
    assert(accel  != 0);

    memset((char*)accel + numPrims * sizeAccelInBytes,0,additional_size); // clear out as a 4-wide access is possible


    bvh4hair->accel = accel;
    bvh4hair->size_node  = size_node;
    bvh4hair->size_accel = size_accel;

    size_prims = size_primrefs;    
  }


  void BVH4HairBuilder::parallelBinningGlobal(const size_t threadID, const size_t numThreads)
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

    Bezier1i  *__restrict__ const tmp_prims = (Bezier1i*)accel;

    fastbin_copy<Bezier1i,false>(prims,tmp_prims,startID,endID,centroidBoundsMin_2,scale,global_bin16[threadID]);    

    scene->lockstep_scheduler.syncThreadsWithReduction( threadID, numThreads, reduceBinsParallel, global_bin16 );
    
    if (threadID == 0)
      {
	const float voxelArea = area(current.bounds.geometry);

	global_sharedData.split.cost = items * voxelArea * INTERSECTION_COST;;
	
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


  void BVH4HairBuilder::parallelPartitioning(BuildRecord& current,
					  Bezier1i * __restrict__ l_source,
					  Bezier1i * __restrict__ r_source,
					  Bezier1i * __restrict__ l_dest,
					  Bezier1i * __restrict__ r_dest,
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

    const mic_m dim_mask = mic_m::shift1[bestSplitDim];

    for (;l_source<r_source;)
      {
	evictL1(l_source-2);	
	    
	prefetch<PFHINT_NT>(l_source+2);
	prefetch<PFHINT_L2>(l_source + L2_PREFETCH_ITEMS + 4);

	const mic2f b = l_source->getBounds();
	const mic_f b_min = b.x;
	const mic_f b_max = b.y;
	// const mic_f b_min = broadcast4to16f(&l_source->lower);
	// const mic_f b_max = broadcast4to16f(&l_source->upper);
	const mic_f b_centroid2 = b_min + b_max;

	if (likely(lt_split(b_min,b_max,dim_mask,c,s,mic_f(bestSplit)))) 
	  {
	    *l_dest++ = *l_source++;

	    leftSceneBoundsMin = min(leftSceneBoundsMin,b_min);
	    leftSceneBoundsMax = max(leftSceneBoundsMax,b_max);

	    leftCentroidBoundsMin = min(leftCentroidBoundsMin,b_centroid2);
	    leftCentroidBoundsMax = max(leftCentroidBoundsMax,b_centroid2);
	    
	  }
	else
	  {
	    *r_dest++ = *l_source++;

	    rightSceneBoundsMin = min(rightSceneBoundsMin,b_min);
	    rightSceneBoundsMax = max(rightSceneBoundsMax,b_max);
	    rightCentroidBoundsMin = min(rightCentroidBoundsMin,b_centroid2);
	    rightCentroidBoundsMax = max(rightCentroidBoundsMax,b_centroid2);
	  }
      }
    
    store4f(&local_left.geometry.lower,leftSceneBoundsMin);
    store4f(&local_left.geometry.upper,leftSceneBoundsMax);
    store4f(&local_left.centroid2.lower,leftCentroidBoundsMin);
    store4f(&local_left.centroid2.upper,leftCentroidBoundsMax);

    store4f(&local_right.geometry.lower,rightSceneBoundsMin);
    store4f(&local_right.geometry.upper,rightSceneBoundsMax);
    store4f(&local_right.centroid2.lower,rightCentroidBoundsMin);
    store4f(&local_right.centroid2.upper,rightCentroidBoundsMax);


  }

  void BVH4HairBuilder::parallelPartitioningGlobal(const size_t threadID, const size_t numThreads)
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

    Bezier1i  *__restrict__ const tmp_prims = (Bezier1i*)accel;

    Bezier1i * __restrict__ l_source = tmp_prims + startID;
    Bezier1i * __restrict__ r_source = tmp_prims + endID;

    Bezier1i * __restrict__ l_dest     = prims + current.begin + thread_start_left;
    Bezier1i * __restrict__ r_dest     = prims + current.begin + thread_start_right + bestSplitNumLeft;

    __aligned(64) Centroid_Scene_AABB local_left;
    __aligned(64) Centroid_Scene_AABB local_right; // just one local

    parallelPartitioning(current,l_source,r_source,l_dest,r_dest,global_sharedData.split,local_left,local_right);

    global_sharedData.left.extend_atomic(local_left); 
    global_sharedData.right.extend_atomic(local_right);  
  }


  void BVH4HairBuilder::parallelBinningLocal(const size_t localThreadID,const size_t globalThreadID)
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

    Bezier1i  *__restrict__ const tmp_prims = (Bezier1i*)accel;

    fastbin_copy<Bezier1i,false>(prims,tmp_prims,startID,endID,centroidBoundsMin_2,scale,global_bin16[globalThreadID]);    

    localTaskScheduler[globalCoreID].syncThreads(localThreadID);

    if (localThreadID == 0)
      {
	Bin16 &bin16 = global_bin16[globalThreadID];

	for (size_t i=1;i<4;i++)
	  bin16.merge(global_bin16[globalThreadID+i]);

	const float voxelArea = area(current.bounds.geometry);

	local_sharedData[globalCoreID].split.cost = items * voxelArea * INTERSECTION_COST;	

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


  void BVH4HairBuilder::parallelPartitioningLocal(const size_t localThreadID,const size_t globalThreadID)
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

    Bezier1i  *__restrict__ const tmp_prims = (Bezier1i*)accel;

    Bezier1i * __restrict__ l_source = tmp_prims + startID;
    Bezier1i * __restrict__ r_source = tmp_prims + endID;

    Bezier1i * __restrict__ l_dest     = prims + current.begin + thread_start_left;
    Bezier1i * __restrict__ r_dest     = prims + current.begin + thread_start_right + bestSplitNumLeft;

    __aligned(64) Centroid_Scene_AABB local_left;
    __aligned(64) Centroid_Scene_AABB local_right; 

    parallelPartitioning(current,l_source,r_source,l_dest,r_dest,sd.split,local_left,local_right);

    sd.left.extend_atomic(local_left); 
    sd.right.extend_atomic(local_right);  
    
  }


  void BVH4HairBuilder::build(const size_t threadIndex, const size_t threadCount) 
  {
    DBG(PING);
    if (threadIndex != 0) {
      FATAL("threadIndex != 0");
    }

    const size_t totalNumPrimitives = getNumPrimitives();


    DBG(DBG_PRINT(totalNumPrimitives));

    /* print builder name */
    if (unlikely(g_verbose >= 2)) 
      printBuilderName();

    if (likely(totalNumPrimitives == 0))
      {
	DBG(std::cout << "EMPTY SCENE BUILD" << std::endl);
	bvh4hair->root = BVH4Hair::invalidNode;
	bvh4hair->bounds = empty;
	bvh4hair->accel = NULL;
	return;
      }


    /* allocate BVH data */
    allocateData(threadCount,totalNumPrimitives);

    if (likely(numPrimitives > SINGLE_THREADED_BUILD_THRESHOLD && threadCount > 1) )
      {
	DBG(std::cout << "PARALLEL BUILD" << std::endl);
	build_main(threadIndex,threadCount);

      }
    else
      {
	/* number of primitives is small, just use single threaded mode */
	assert( numPrimitives > 0 );
	DBG(std::cout << "SERIAL BUILD" << std::endl);
	build_main(0,1);
      }

    if (g_verbose >= 2) {
      double perf = totalNumPrimitives/dt*1E-6;
      std::cout << "[DONE] " << 1000.0f*dt << "ms (" << perf << " Mtris/s), primitives " << numPrimitives << std::endl;
      std::cout << BVH4HairStatistics<BVH4Hair::UnalignedNode>(bvh4hair).str();
    }

  }




  void BVH4HairBuilder::build_main(size_t threadIndex, size_t threadCount) 
  {
    DBG(PING);

    TIMER(double msec = 0.0);

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
    atomicID.reset(1);
    node[0].setInvalid();
    node[0].setMatrix(global_bounds.geometry,0);

    /* create initial build record */
    BuildRecord br;
    br.init(global_bounds,0,numPrimitives);
    br.depth       = 1;
    //br.parentID    = 0;
    br.parentPtr   = &node[0].child(0);

    /* node allocator */
    NodeAllocator alloc(atomicID,numAllocated64BytesBlocks);
        
    /* push initial build record to global work stack */
    global_workStack.reset();
    global_workStack.push_nolock(br);    

    /* work in multithreaded toplevel mode until sufficient subtasks got generated */    
    const size_t coreCount = (threadCount+3)/4;
    while (global_workStack.size() < coreCount &&
	   global_workStack.size()+BVH4Hair::N <= SIZE_GLOBAL_WORK_STACK) 
    {
      BuildRecord br;
      if (!global_workStack.pop_nolock_largest(br)) break;
      recurseSAH(br,alloc,BUILD_TOP_LEVEL,threadIndex,threadCount);      
    }

    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "build_top_level " << 1000. * msec << " ms" << std::endl << std::flush);

    /* fill per core work queues */    
    TIMER(msec = getSeconds());    
    scene->lockstep_scheduler.dispatchTask(task_fillLocalWorkQueues, this, threadIndex, threadCount );
    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "task_fillLocalWorkQueues " << 1000. * msec << " ms" << std::endl << std::flush);

    /* now process all created subtasks on multiple threads */    
    TIMER(msec = getSeconds());    
    scene->lockstep_scheduler.dispatchTask(task_buildSubTrees, this, threadIndex, threadCount );
    DBG(DBG_PRINT(atomicID));
    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "task_buildSubTrees " << 1000. * msec << " ms" << std::endl << std::flush);

    numNodes = atomicID; 
    
    /* update BVH4 */
    
    bvh4hair->root             = node[0].child(0);
    bvh4hair->bounds           = global_bounds.geometry;
    bvh4hair->unaligned_nodes  = (BVH4Hair::UnalignedNode*)node;
    bvh4hair->accel            = prims;

    /* stop measurement */
    if (g_verbose >= 2) 
      dt = getSeconds()-t0;
  }

  bool split_fallback(Bezier1i * __restrict__ const primref, BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild)
  {
    // unsigned int blocks4 = (current.items()+3)/4;
    // unsigned int center = current.begin + (blocks4/2)*4; 

    unsigned int center = (current.begin + current.end)/2;
    assert(center != current.begin);
    assert(center != current.end);

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

  __forceinline bool BVH4HairBuilder::split(BuildRecord& current, BuildRecord& left, BuildRecord& right, const size_t mode, const size_t threadID, const size_t numThreads)
  {
    DBG(PING);

   if (unlikely(mode == BUILD_TOP_LEVEL))
      {
	DBG(DBG_PRINT("TOP_LEVEL"));

	if (current.items() >= BUILD_RECORD_PARALLEL_SPLIT_THRESHOLD && numThreads > 1)
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


  __forceinline bool BVH4HairBuilder::splitSequential(BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild)
  {
    /* mark as leaf if leaf threshold reached */
    if (current.items() <= MAX_ITEMS_PER_LEAF) {
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

    fastbin<Bezier1i>(prims,current.begin,current.end,centroidBoundsMin_2,scale,leftArea,rightArea,leftNum);

    const unsigned int items = current.items();
    const float voxelArea = area(current.bounds.geometry);
    Split split;
    split.cost = items * voxelArea * INTERSECTION_COST;;

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

	const unsigned int mid = partitionPrimitives<L2_PREFETCH_ITEMS>(prims ,current.begin, current.end-1, split.pos, split.dim, centroidBoundsMin_2, scale, leftChild.bounds, rightChild.bounds);

	assert(area(leftChild.bounds.geometry) >= 0.0f);
	//assert(current.begin + mid == current.begin + split.numLeft) // can happen

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



    if (leftChild.items()  <= MAX_ITEMS_PER_LEAF) leftChild.createLeaf();
    if (rightChild.items() <= MAX_ITEMS_PER_LEAF) rightChild.createLeaf();	
    return true;
  }
  

  bool BVH4HairBuilder::splitParallelGlobal( BuildRecord &current,
					  BuildRecord &leftChild,
					  BuildRecord &rightChild,
					  const size_t threadID,
					  const size_t numThreads)
  {
    DBG(PING);

    const unsigned int items = current.end - current.begin;
    assert(items >= BUILD_RECORD_PARALLEL_SPLIT_THRESHOLD);
  
    /* mark as leaf if leaf threshold reached */
    if (items <= MAX_ITEMS_PER_LEAF) {
      current.createLeaf();
      return false;
    }

     global_sharedData.rec = current;
     global_sharedData.split.reset();
     global_sharedData.left.reset();
     global_sharedData.right.reset();
     
     scene->lockstep_scheduler.dispatchTask( task_parallelBinningGlobal, this, threadID, numThreads );

     if (unlikely(global_sharedData.split.pos == -1)) 
       split_fallback(prims,current,leftChild,rightChild);
     else
       {
	 global_sharedData.left.reset();
	 global_sharedData.right.reset();

	 global_sharedData.lCounter.reset(0);
	 global_sharedData.rCounter.reset(0); 

	 scene->lockstep_scheduler.dispatchTask( task_parallelPartitioningGlobal, this, threadID, numThreads );

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
     
     if (leftChild.items()  <= MAX_ITEMS_PER_LEAF) leftChild.createLeaf();
     if (rightChild.items() <= MAX_ITEMS_PER_LEAF) rightChild.createLeaf();
     return true;
  }


  bool BVH4HairBuilder::splitParallelLocal(BuildRecord &current,
					   BuildRecord &leftChild,
					   BuildRecord &rightChild,
					   const size_t threadID)
  {
    const unsigned int items    = current.end - current.begin;
    const size_t globalCoreID   = threadID / 4;
    const size_t localThreadID  = threadID % 4;
    const size_t globalThreadID = threadID;
    
    assert(items >= THRESHOLD_FOR_SUBTREE_RECURSION);
  
    /* mark as leaf if leaf threshold reached */
    if (items <= MAX_ITEMS_PER_LEAF) {
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
     
     if (leftChild.items()  <= MAX_ITEMS_PER_LEAF) leftChild.createLeaf();
     if (rightChild.items() <= MAX_ITEMS_PER_LEAF) rightChild.createLeaf();
     return true;
  }


  __forceinline void BVH4HairBuilder::createLeaf(BuildRecord& current, NodeAllocator& alloc,const size_t threadIndex, const size_t threadCount)
  {
#if defined(DEBUG)
    if (current.depth > BVH4Hair::maxBuildDepthLeaf) 
      THROW_RUNTIME_ERROR("ERROR: depth limit reached");
#endif
    
    /* create leaf */
    if (current.items() <= MAX_ITEMS_PER_LEAF) {
      createBVH4HairLeaf(current.parentPtr,current.begin,current.items());
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
    size_t numChildren = 1;
    const size_t currentIndex = alloc.get(1);

    createBVH4HairNode(current.parentPtr,currentIndex);
    
    node[currentIndex].prefetchNode<PFHINT_L2EX>();

    /* recurse into each child */
    for (size_t i=0; i<numChildren; i++) 
    {
      node[currentIndex].setMatrix(children[i].bounds.geometry, i);
      children[i].parentPtr = &node[currentIndex].child(i);
      children[i].depth     = current.depth+1;
      createLeaf(children[i],alloc,threadIndex,threadCount);
    }
  }  


  void BVH4HairBuilder::recurseSAH(BuildRecord& current, 
				   NodeAllocator& alloc,
				   const size_t mode, 
				   const size_t threadID, 
				   const size_t numThreads)
  {
    __aligned(64) BuildRecord children[BVH4Hair::N];

    /* create leaf node */
    if (current.depth >= BVH4Hair::maxBuildDepth || current.isLeaf()) {
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
      
    } while (numChildren < BVH4Hair::N);

    /* create leaf node if no split is possible */
    if (numChildren == 1) {
      createLeaf(current,alloc,threadID,numThreads);
      return;
    }

    /* allocate next four nodes */
    const size_t currentIndex = alloc.get(1);
    
    BVH4Hair::UnalignedNode *current_node = (BVH4Hair::UnalignedNode *)&node[currentIndex];

#if ENABLE_AABB_NODES == 1
    createBVH4HairNode(current.parentPtr,currentIndex,BVH4Hair::alignednode_mask);
#else
    createBVH4HairNode(current.parentPtr,currentIndex);
#endif


    /* init used/unused nodes */
    current_node->prefetchNode<PFHINT_L2EX>();
    current_node->setInvalid();

    /* recurse into each child */
    for (unsigned int i=0; i<numChildren; i++) 
    {
      current_node->setMatrix(children[i].bounds.geometry,i);
      //children[i].parentID    = currentIndex;
      children[i].parentPtr   = &current_node->child(i);      
      recurse(children[i],alloc,mode,threadID,numThreads);
    }    

    DBG(DBG_PRINT(*current_node));

  }


  __forceinline void BVH4HairBuilder::recurse(BuildRecord& current, NodeAllocator& alloc,const size_t mode, const size_t threadID, const size_t numThreads)
  {
    if (mode == BUILD_TOP_LEVEL) {
      global_workStack.push_nolock(current);
    }
    else if (mode == FILL_LOCAL_QUEUES) {
      const size_t coreID = threadID/4;
      if (!local_workStack[coreID].push(current))
        recurseSAH(current,alloc,RECURSE,threadID,numThreads);
    }
    else
      {
	BuildRecordOBB current_obb;
	current_obb = current;

	computeUnalignedSpace(current_obb);
	computeUnalignedSpaceBounds(current_obb);
   
#if defined(THRESHOLD_SWITCH)
	if (area( current_obb.bounds.geometry ) < area( current.bounds.geometry ) * AABB_OBB_SWITCH_THRESHOLD)
	  recurseOBB(current_obb,alloc,/*mode*/ RECURSE,threadID,numThreads);
	else
	  recurseSAH(current,alloc,RECURSE,threadID,numThreads);
#else
	recurseOBB(current_obb,alloc,/*mode*/ RECURSE,threadID,numThreads);
#endif


      }
  }

  // ===============================================================================================================================================
  // ===============================================================================================================================================
  // ===============================================================================================================================================


  void BVH4HairBuilder::recurseOBB(BuildRecordOBB& current, NodeAllocator& alloc, const size_t mode, const size_t threadID, const size_t numThreads)
  {
#if ENABLE_OBB_BVH4 == 0
    FATAL("recurseOBB disabled");
#endif

    __aligned(64) BuildRecordOBB children[BVH4Hair::N];

    /* create leaf node */
    if (current.depth >= BVH4Hair::maxBuildDepth || current.isLeaf()) {
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
      __aligned(64) BuildRecordOBB left, right;
      if (!splitSequentialOBB(children[bestChild],left,right,false)) 
        continue;
      
      /* add new children left and right */
      left.depth = right.depth = current.depth+1;
      children[bestChild] = children[numChildren-1];
      children[numChildren-1] = left;
      children[numChildren+0] = right;
      numChildren++;
      
    } while (numChildren < BVH4Hair::N);

    /* create leaf node if no split is possible */
    if (numChildren == 1) {
      createLeaf(current,alloc,threadID,numThreads);
      return;
    }

    /* allocate next four nodes */
    const size_t currentIndex = alloc.get(1);
    /* recurseOBB */

    node[currentIndex].prefetchNode<PFHINT_L2EX>();

    /* init used/unused nodes */
    node[currentIndex].setInvalid();

    // === default OBB node ===
    createBVH4HairNode(current.parentPtr,currentIndex);

    for (unsigned int i=0; i<numChildren; i++) 
      node[currentIndex].setMatrix(children[i].xfm,children[i].bounds.geometry,i);

    /* recurse into each child */
    for (unsigned int i=0; i<numChildren; i++) 
      {
	children[i].parentPtr = &node[currentIndex].child(i);
	recurseOBB(children[i],alloc,mode,threadID,numThreads);
      }  

  }

  void BVH4HairBuilder::buildSubTree(BuildRecord& current, 
				     NodeAllocator& alloc, 
				     const size_t mode,
				     const size_t threadID, 
				     const size_t numThreads)
  {
    DBG(PING);
#if ENABLE_OBB_BVH4 == 1

    if (mode == FILL_LOCAL_QUEUES)
      {
	recurseSAH(current,alloc,mode,threadID,numThreads);
	return;
      }
    else
      {

	BuildRecordOBB current_obb;
	current_obb = current;

	computeUnalignedSpace(current_obb);
	computeUnalignedSpaceBounds(current_obb);

#if defined( THRESHOLD_SWITCH ) 
	if (area( current_obb.bounds.geometry ) < area( current.bounds.geometry ) * AABB_OBB_SWITCH_THRESHOLD)
	  recurseOBB(current_obb,alloc,/*mode*/ RECURSE,threadID,numThreads);
	else
	  recurseSAH(current,alloc,/*mode*/ RECURSE,threadID,numThreads);
#else
	recurseOBB(current_obb,alloc,/*mode*/ RECURSE,threadID,numThreads);

#endif
      }

#else
    recurseSAH(current,alloc,/*mode*/ RECURSE,threadID,numThreads);
#endif
  }


  __forceinline bool BVH4HairBuilder::splitSequentialOBB(BuildRecordOBB& current, BuildRecordOBB& leftChild, BuildRecordOBB& rightChild, const bool binAABB)
  {
    DBG(PING);
    DBG(DBG_PRINT(current));

    //computeUnalignedSpace(current);
    //computeUnalignedSpaceBounds(current);

    /* mark as leaf if leaf threshold reached */
    if (current.items() <= MAX_ITEMS_PER_LEAF) {
      current.createLeaf();
      return false;
    }


    const unsigned int items = current.items();
    const float voxelArea = area(current.bounds.geometry); 

    Split split;
    split.cost = items * voxelArea * INTERSECTION_COST;;
    
    const mic_f centroidMin = broadcast4to16f(&current.bounds.centroid2.lower);
    const mic_f centroidMax = broadcast4to16f(&current.bounds.centroid2.upper);

    const mic_f centroidBoundsMin_2 = centroidMin;
    const mic_f centroidDiagonal_2  = centroidMax-centroidMin;
    const mic_f scale = select(centroidDiagonal_2 != 0.0f,rcp(centroidDiagonal_2) * mic_f(16.0f * 0.99f),mic_f::zero());

    mic_f leftArea[3];
    mic_f rightArea[3];
    mic_i leftNum[3];


    const mic3f cmat = convert(current.xfm);
    
    fastbin_xfm<Bezier1i>(prims,cmat,current.begin,current.end,centroidBoundsMin_2,scale,leftArea,rightArea,leftNum);


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

    if (binAABB)
      {
	Split splitAABB;
	splitAABB.cost = items * voxelArea * INTERSECTION_COST;

	fastbin<Bezier1i>(prims,current.begin,current.end,centroidBoundsMin_2,scale,leftArea,rightArea,leftNum);

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

	    if (lt(cost,mic_f(splitAABB.cost)))
	      {
		const mic_f min_cost    = vreduce_min(cost); 
		const mic_m m_pos       = min_cost == cost;
		const unsigned long pos = bitscan64(m_pos);	    
		assert(pos < 15);
		if (pos < 15)
		  {
		    splitAABB.cost    = cost[pos];
		    splitAABB.pos     = pos+1;
		    splitAABB.dim     = dim;	    
		    splitAABB.numLeft = lnum[pos];
		  }
	      }
	  };

      }

    if (unlikely(split.pos == -1)) 
      {
	split_fallback(prims,current,leftChild,rightChild);
      }
    else 
      {
	/* partitioning of items */
	leftChild.bounds.reset();
	rightChild.bounds.reset();

	const unsigned int mid = partitionPrimitives_xfm<L2_PREFETCH_ITEMS>(prims,cmat,current.begin, current.end-1, split.pos, split.dim, centroidBoundsMin_2, scale, leftChild.bounds, rightChild.bounds);

	assert(area(leftChild.bounds.geometry) >= 0.0f);
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

    computeUnalignedSpace(leftChild);
    computeUnalignedSpaceBounds(leftChild);

    computeUnalignedSpace(rightChild);
    computeUnalignedSpaceBounds(rightChild);

    //leftChild.xfm  = current.xfm;
    //rightChild.xfm = current.xfm;

    //if (leftChild.items()  <= MAX_ITEMS_PER_LEAF) leftChild.createLeaf();
    //if (rightChild.items() <= MAX_ITEMS_PER_LEAF) rightChild.createLeaf();

    return true;
  }


  __forceinline void BVH4HairBuilder::computeUnalignedSpace( BuildRecordOBB& current )
  {
    Vec3fa axis(0,0,1);
#if 1
    for (size_t i=current.begin;i<current.end;i++)
      {
	const Bezier1i &b = prims[i];
	const Vec3fa &p0 = b.p[0];
	const Vec3fa &p3 = b.p[3];
	const Vec3fa axis1 = normalize(p3 - p0);
	if (length(p3 - p0) > 1E-9f) {
	  axis = axis1;
	  break;
	}	
      }
#endif
    current.xfm = frame(axis).transposed();    

    current.PreQuantizeMatrix();

#if 0
    PING;
    DBG_PRINT( current.xfm );
    DBG_PRINT( current.xfm );

#endif
    DBG(DBG_PRINT(current.xfm));
  }

  __forceinline void BVH4HairBuilder::computeUnalignedSpaceBounds( BuildRecordOBB& current )
  {
    const mic_f c0 = broadcast4to16f((float*)&current.xfm.vx);
    const mic_f c1 = broadcast4to16f((float*)&current.xfm.vy);
    const mic_f c2 = broadcast4to16f((float*)&current.xfm.vz);

    const mic_f p_inf( pos_inf );
    const mic_f n_inf( neg_inf );

    mic2f centroid2(p_inf, n_inf);
    mic2f geometry(p_inf, n_inf);

    for (size_t i=current.begin;i<current.end;i++)
      {
	prefetch<PFHINT_NT>(&prims[i+4]);
	prefetch<PFHINT_L2>(&prims[i+16]);
	const mic2f b = prims[i].getBounds(c0,c1,c2);
	const mic_f b_min = b.x;
	const mic_f b_max = b.y;
	const mic_f c2    = b_min + b_max;
	centroid2.x = min(centroid2.x, c2);
	centroid2.y = max(centroid2.y, c2);
	geometry.x  = min(geometry.x, b_min); 
	geometry.y  = max(geometry.y, b_max); 	
      }    
    
    store4f(&current.bounds.centroid2.lower,centroid2.x);
    store4f(&current.bounds.centroid2.upper,centroid2.y);

    store4f(&current.bounds.geometry.lower,geometry.x);
    store4f(&current.bounds.geometry.upper,geometry.y);

    DBG(
	for (size_t i=current.begin;i<current.end;i++)
	  {
	    for (size_t j=0;j<4;j++)
	      {
		const Vec3fa v = prims[i].p[j];
		const Vec3fa xfm_v = xfmPoint(current.xfm,v);
		if (disjoint(current.bounds.geometry,xfm_v))
		  {
		    DBG_PRINT(v);
		    DBG_PRINT(xfm_v);
		    DBG_PRINT(current.xfm);
		    DBG_PRINT(current.bounds.geometry);
		    exit(0);
		  }
	      }
	  }
	
	);
    current.sArea = area(current.bounds.geometry);

  }


  void BVH4HairBuilder::createAccel(const size_t threadIndex, const size_t threadCount)
  {
    DBG(PING);
  }

  std::string BVH4HairBuilder::getStatistics()
  {
    return std::string("not implemented");
  }


};
