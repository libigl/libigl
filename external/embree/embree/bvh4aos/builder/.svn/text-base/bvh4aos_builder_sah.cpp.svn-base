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

#include "common/default.h"
#include "bvh4aos_globals.h"
#include "bvh4aos_builder_util.h"
#include "bvh4aos_builder_common.h"

#define TIMER(x) 
//#define CHECK(x) { LOCK; x ; UNLOCK; }
#define CHECK(x) 

#define PRE_ALLOCATION_FACTOR 0.6f
#define L2_PREFETCH_DISTANCE  16

#define LOCAL_BUILD_RECORD_STACK_PUSH_THRESHOLD 128
#define BUILD_RECORD_PARALLEL_SPLIT_THRESHOLD   256
#define MAX_BUILD_RECORD_STACK_ENTRIES          ParallelQBVHBuilder::numActiveLocalWorkQueues * 2
#define BUILD_RECORD_SPLIT_THRESHOLD            ParallelQBVHBuilder::numActiveLocalWorkQueues * 16

namespace embree
{

  // =================================================================
  // =================================================================
  // =================================================================


  AABB            * __restrict__ ParallelQBVHBuilder::aabb                   = NULL; 
  AABB            * __restrict__ ParallelQBVHBuilder::tmp_aabb               = NULL;
  BVHNode         * __restrict__ ParallelQBVHBuilder::globalNodePtr          = NULL;
  Vec3fa          * __restrict__ ParallelQBVHBuilder::globalVertexPtr        = NULL;  
  BuilderTriangle * __restrict__ ParallelQBVHBuilder::globalBuildTrianglePtr = NULL;
  TriangleAccel   * __restrict__ ParallelQBVHBuilder::globalAccelPtr         = NULL;

  unsigned int (* ParallelQBVHBuilder::recursePtr)(BuildRecord &current, 
						   const unsigned int mode, 
						   const unsigned int threadID,
						   unsigned int &localNodeID,
						   unsigned int &localNodeIDs) = NULL;

  size_t ParallelQBVHBuilder::triangles  = 0;
  size_t ParallelQBVHBuilder::vertices   = 0;
  size_t ParallelQBVHBuilder::qbvh_nodes = 0;
  size_t ParallelQBVHBuilder::max_nodes  = 0;
  unsigned int ParallelQBVHBuilder::numActiveLocalWorkQueues = 0;

  bool ParallelQBVHBuilder::enableTaskStealing = true;

  MIC_ALIGN AtomicCounter  ParallelQBVHBuilder::atomicID  = 0;
  MIC_ALIGN AABB ParallelQBVHBuilder::initQBVHNode[2];
  MIC_ALIGN Centroid_Scene_AABB ParallelQBVHBuilder::global_cs_AABB;
  MIC_ALIGN WorkStack<BuildRecord,NUM_ATOMIC_BUILD_RECORDS> ParallelQBVHBuilder::atomicBuildRecordStack;
  MIC_ALIGN WorkStack<BuildRecord,NUM_LOCAL_ATOMIC_BUILD_RECORDS> ParallelQBVHBuilder::localAtomicBuildRecordStack[MAX_MIC_CORES];
  MIC_ALIGN AtomicMutex ParallelQBVHBuilder::threadMutex;
  MIC_ALIGN ParallelQBVHBuilder::SplitData ParallelQBVHBuilder::splitData;
  MIC_ALIGN Bin16 ParallelQBVHBuilder::threadBin16[MAX_MIC_THREADS];



  static MIC_ALIGN AtomicCounter lCounter;
  static MIC_ALIGN AtomicCounter rCounter;

  // =================================================================
  // =================================================================
  // =================================================================




  // =================================================================
  // =================================================================
  // =================================================================


  void ParallelQBVHBuilder::thread_build_local(const unsigned int threadID)
  {
    MIC_ALIGN BuildRecord br;

    const unsigned int coreID = threadID >> 2;
    const unsigned int localThreadID = threadID % 4;

    localAtomicBuildRecordStack[coreID].prefetchL1Ex();
    localAtomicBuildRecordStack[coreID].get_addr(threadID%4)->prefetchL1Ex();

    unsigned int localNodeID = ParallelQBVHBuilder::atomicID.add(LOCAL_NODE_IDS);
    unsigned int localNodeIDs = 0;

    CHECK(
	  std::cout << "threadID " << threadID << " coreID " << coreID << " slots " << localAtomicBuildRecordStack[coreID].usedSlots() << std::endl << std::flush
	  );


    // =======================================================================================

    // stay local on the core as long as possible to maximize cache reusage 

    while(1)
      {
	if ( localAtomicBuildRecordStack[coreID].pop_smallest(br,threadID) == true)
	  {
	    recursePtr(br,PUSH_LOCAL,threadID,localNodeID,localNodeIDs);
	  }
	else
	  {
	    delayThread(1024);
	    if (unlikely(localAtomicBuildRecordStack[coreID].allThreadsDone())) break;
	  }
      }

    // =======================================================================================

    while(1)
      {
	while(1)
	  {
	    if ( localAtomicBuildRecordStack[coreID].pop_smallest(br,threadID) == true)
	      {
		recursePtr(br,PUSH_LOCAL,threadID,localNodeID,localNodeIDs);
	      }
	    else
	      {
		break;
	      }
	  }
	if (enableTaskStealing)
	  {
	    // === task stealing ===
	    unsigned int stillWorking = 0;
	    int index = -1;
	    for (unsigned int i=0;i<numActiveLocalWorkQueues;i++)
	      {
		unsigned int newCoreID = coreID + i;
		if (newCoreID >= numActiveLocalWorkQueues) newCoreID -= numActiveLocalWorkQueues;

		stillWorking |= localAtomicBuildRecordStack[newCoreID].getAllThreadStates();
		if (localAtomicBuildRecordStack[newCoreID].usedSlots())
		  if (localAtomicBuildRecordStack[newCoreID].pop(br) == true)
		    {
		      index = newCoreID;
		      break;
		    }
	      }
	    if (index >= 0)
	      {
		localAtomicBuildRecordStack[coreID].setThreadState(threadID,1); 
		recursePtr(br,PUSH_LOCAL,threadID,localNodeID,localNodeIDs);
		localAtomicBuildRecordStack[coreID].setThreadState(threadID,0); 
	      }
	    else
	      {
		if (stillWorking == 0)
		  break;
		else
		  delayThread(1024);
	      }
	  }
	else
	  break;
	// =====================
      }

  }

  // =================================================================
  // =================================================================
  // =================================================================


  void ParallelQBVHBuilder::thread_createAABBs(const unsigned int threadID)
  {
    if (unlikely(threadID >= splitData.worker_threads)) return;

    const unsigned int WORKERS  = splitData.worker_threads;
    const unsigned int workerID = threadID;

    unsigned int primitives_per_block = triangles / WORKERS;
    if (primitives_per_block % 2) primitives_per_block -= primitives_per_block % 2;
    const unsigned int startID = workerID * primitives_per_block;
    const unsigned int endID  = (workerID != WORKERS-1) ? (startID + primitives_per_block) : triangles - (triangles % 2);


    if (unlikely(startID >= endID)) return;

    assert(startID % 2 == 0);
    assert(endID   % 2 == 0);
    assert(endID <= triangles);

    MIC_ALIGN Centroid_Scene_AABB cs_AABB;
    cs_AABB.reset();

    const Vec3fa          * __restrict__ const vptr = globalVertexPtr;
    const BuilderTriangle * __restrict__       tptr = globalBuildTrianglePtr + startID;
    AABB                  * __restrict__       aptr = aabb + startID;

    prefetch<PFHINT_L1>(tptr);


    mic_i triID = mic_i(startID);
    int tt = startID;
    for (unsigned int i=startID;i<endID;i+=2,tptr+=2,aptr+=2, tt+=2)
      {
	evictL2(tptr-4);

	/* --- assume item is identity --- */

	prefetch<PFHINT_L2>(&vptr[tptr[1].v0]);
	prefetch<PFHINT_L2>(&vptr[tptr[1].v1]);
	prefetch<PFHINT_L2>(&vptr[tptr[1].v2]);

	prefetch<PFHINT_L2>(&vptr[tptr[0].v1]);
	prefetch<PFHINT_L2>(&vptr[tptr[0].v2]);

	const mic_f v0_0 = upconv4f((float*)&vptr[tptr[0].v0]); 
	const mic_f v1_0 = upconv4f((float*)&vptr[tptr[0].v1]);
	const mic_f v2_0 = upconv4f((float*)&vptr[tptr[0].v2]);

	const mic_f b_min_0 = _min(v0_0,_min(v1_0,v2_0));
	const mic_f b_max_0 = _max(v0_0,_max(v1_0,v2_0));

	const mic_f centroid_0 = (b_min_0+b_max_0) * mic_f(0.5f); 

	cs_AABB.extend(b_min_0,b_max_0,centroid_0);
      
	const mic_f b0 = sel(0x8888,cast_to_mic_f(triID),b_min_0);
	const mic_f b1 = b_max_0;

	triID += mic_i::one();

	prefetch<PFHINT_L1>(tptr + 1);

	const mic_f v0_1 = upconv4f((float*)&vptr[tptr[1].v0]); 
	const mic_f v1_1 = upconv4f((float*)&vptr[tptr[1].v1]);
	const mic_f v2_1 = upconv4f((float*)&vptr[tptr[1].v2]);

	const mic_f b_min_1 = _min(v0_1,_min(v1_1,v2_1));
	const mic_f b_max_1 = _max(v0_1,_max(v1_1,v2_1));

	const mic_f centroid_1 = (b_min_1+b_max_1) * mic_f(0.5f); 

	cs_AABB.extend(b_min_1,b_max_1,centroid_1);
      

	const mic_f b2 = sel(0x8888,cast_to_mic_f(triID),b_min_1);
	const mic_f b3 = b_max_1;

	triID += mic_i::one();

	const mic_f twoAABBs = lane_shuffle_gather<0>(b0,b1,b2,b3);

	store16f_ngo(aptr[0].m_min,twoAABBs);                              
      }  

    prefetch<PFHINT_L1EX>(&global_cs_AABB);

    if (workerID == WORKERS-1)
      if (triangles%2)
	{
	  prefetch<PFHINT_L1EX>(aptr);
	  const mic_f v0 = upconv4f((float*)&vptr[tptr->v0]); 
	  const mic_f v1 = upconv4f((float*)&vptr[tptr->v1]);
	  const mic_f v2 = upconv4f((float*)&vptr[tptr->v2]);
	  const mic_f b_min = _min(v0,_min(v1,v2));
	  const mic_f b_max = _max(v0,_max(v1,v2));
	  const mic_f centroid = (b_min+b_max) * mic_f(0.5f); // (v0 + v1 + v2) * const1_3;
	  store4f(aptr->m_max,b_max); 
	  store4f(aptr->m_min,b_min); 
	  cs_AABB.extend(b_min,b_max,centroid);
	  aptr->ext_min.t = triangles-1;
	  evictL1(aptr);
	}  

    global_cs_AABB.extend_atomic(cs_AABB);
  }

  // =================================================================
  // =================================================================
  // =================================================================

  template<unsigned int DISTANCE>
  _INLINE unsigned int partitionAABBs(AABB *__restrict__ aabb,
				      const unsigned int begin,
				      const unsigned int end,
				      const unsigned int bestSplit,
				      const unsigned int bestSplitDim,
				      const mic_f &centroidBoundsMin_2,
				      const mic_f &scale,
				      Centroid_Scene_AABB & local_left,
				      Centroid_Scene_AABB & local_right)
  {
    assert(begin <= end);

    AABB *__restrict__ l = aabb + begin;
    AABB *__restrict__ r = aabb + end;

    const mic_f c = mic_f(centroidBoundsMin_2[bestSplitDim]);
    const mic_f s = mic_f(scale[bestSplitDim]);

    mic_f left_centroidMinAABB = local_left.minCentroidAABB();
    mic_f left_centroidMaxAABB = local_left.maxCentroidAABB();
    mic_f left_sceneMinAABB    = local_left.minSceneAABB();
    mic_f left_sceneMaxAABB    = local_left.maxSceneAABB();

    mic_f right_centroidMinAABB = local_right.minCentroidAABB();
    mic_f right_centroidMaxAABB = local_right.maxCentroidAABB();
    mic_f right_sceneMinAABB    = local_right.minSceneAABB();
    mic_f right_sceneMaxAABB    = local_right.maxSceneAABB();

    const mic_f bestSplit_f = mic_f(bestSplit);
    while(1)
      {
	while (likely(l < r && lt_split(l,bestSplitDim,c,s,bestSplit_f))) 
	  {
	    prefetch<PFHINT_L2EX>(l + DISTANCE);	  
	    {
	      const mic_f b_min = upconv4f((float*)l->m_min);
	      const mic_f b_max = upconv4f((float*)l->m_max);
	      const mic_f centroid = (b_min+b_max) * 0.5f;
	      left_centroidMinAABB = _min(left_centroidMinAABB,centroid);
	      left_centroidMaxAABB = _max(left_centroidMaxAABB,centroid);
	      left_sceneMinAABB    = _min(left_sceneMinAABB,b_min);
	      left_sceneMaxAABB    = _max(left_sceneMaxAABB,b_max);
	    }
	    //evictL1(l-2);

	    ++l;
	  }
	while (likely(l < r && ge_split(r,bestSplitDim,c,s,bestSplit_f))) 
	  {
	    prefetch<PFHINT_L2EX>(r - DISTANCE);
	    {
	      const mic_f b_min = upconv4f((float*)r->m_min);
	      const mic_f b_max = upconv4f((float*)r->m_max);
	      const mic_f centroid = (b_min+b_max) * 0.5f;
	      right_centroidMinAABB = _min(right_centroidMinAABB,centroid);
	      right_centroidMaxAABB = _max(right_centroidMaxAABB,centroid);
	      right_sceneMinAABB    = _min(right_sceneMinAABB,b_min);
	      right_sceneMaxAABB    = _max(right_sceneMaxAABB,b_max);
	    }
	    //evictL1(r+2);


	    --r;
	  }
	if (unlikely(l == r)) {
	  if ( ge_split(r,bestSplitDim,c,s,bestSplit_f))
	    {
	      {
		const mic_f b_min = upconv4f((float*)r->m_min);
		const mic_f b_max = upconv4f((float*)r->m_max);
		const mic_f centroid = (b_min+b_max) * 0.5f;
		right_centroidMinAABB = _min(right_centroidMinAABB,centroid);
		right_centroidMaxAABB = _max(right_centroidMaxAABB,centroid);
		right_sceneMinAABB    = _min(right_sceneMinAABB,b_min);
		right_sceneMaxAABB    = _max(right_sceneMaxAABB,b_max);
	      }	    
	      //local_right.extend(*r);
	    }
	  else 
	    l++; 
	  break;
	}

	xchg(*l,*r);
      }

    local_left.setMinCentroidAABB(left_centroidMinAABB);
    local_left.setMaxCentroidAABB(left_centroidMaxAABB);
    local_left.setMinSceneAABB(left_sceneMinAABB);
    local_left.setMaxSceneAABB(left_sceneMaxAABB);

    local_right.setMinCentroidAABB(right_centroidMinAABB);
    local_right.setMaxCentroidAABB(right_centroidMaxAABB);
    local_right.setMinSceneAABB(right_sceneMinAABB);
    local_right.setMaxSceneAABB(right_sceneMaxAABB);

    assert( aabb + begin <= l && l <= aabb + end);
    assert( aabb + begin <= r && r <= aabb + end);

    return l - (aabb + begin);
  }


  // =================================================================
  // =================================================================
  // =================================================================


  bool ParallelQBVHBuilder::split(BuildRecord &current,
				  BuildRecord &leftChild,
				  BuildRecord &rightChild)
  {
    MIC_ALIGN Centroid_Scene_AABB left;
    MIC_ALIGN Centroid_Scene_AABB right;

    const unsigned int items = current.end - current.begin;
    const float voxelArea    = current.cs_AABB.sceneArea();
    const float centroidArea = current.cs_AABB.centroidArea();
    assert(items > 0);
  
    if (items <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD) 
      {
      createLeaf:      
	current.createLeaf();
	return false;
      }

    if (unlikely(voxelArea <= 0.0f ||
		 eqz(0x7777,current.cs_AABB.centroidDiagonal()) == 0x7777)) // special cases
      {
	assert(items > 4);
	const unsigned int current_mid = current.begin + items/2;
	leftChild.init(current.cs_AABB,current.begin,current_mid);
	rightChild.init(current.cs_AABB,current_mid,current.end);
	return true;
      }

    /* --------------------- */
    /* --- bin centroids --- */
    /* --------------------- */

    const mic_f centroidBoundsMin_2 = current.cs_AABB.minCentroidAABB()  * 2.0f;
    const mic_f centroidDiagonal_2  = current.cs_AABB.centroidDiagonal() * 2.0f;
    const mic_f scale = sel(nz(centroidDiagonal_2),rcp(centroidDiagonal_2) * mic_f(16.0f * 0.99f),mic_f::zero());

    //MIC_ALIGN Bin16 bin16;  

    mic_f leftArea[3];
    mic_f rightArea[3];
    mic_i leftNum[3];

    fastbin(aabb,current.begin,current.end,centroidBoundsMin_2,scale,leftArea,rightArea,leftNum);

    int bestSplit = -1;
    int bestSplitDim = -1;
    int bestNumLeft = -1;
    float bestCost = items * voxelArea;

    for (unsigned int dim = 0;dim < 3;dim++) 
      {
	if (unlikely(centroidDiagonal_2[dim] == 0.0f)) continue;

	const mic_f rArea   = rightArea[dim]; // bin16.prefix_area_rl(dim);
	const mic_f lArea   = leftArea[dim];  // bin16.prefix_area_lr(dim);      
	const mic_i lnum    = leftNum[dim];   // bin16.prefix_count(dim);

	const mic_i rnum    = mic_i(items) - lnum;
	const mic_i lblocks = (lnum + mic_i(3)) >> 2;
	const mic_i rblocks = (rnum + mic_i(3)) >> 2;
	const mic_m m_lnum  = eqz(lnum);
	const mic_m m_rnum  = eqz(rnum);
	const mic_f cost    = sel(m_lnum|m_rnum,mic_f::inf(),lArea * mic_f(lblocks) + rArea * mic_f(rblocks) + voxelArea );
	if (lt(cost,bestCost))
	  {

	    const mic_f min_cost    = set_min16(cost); 
	    const mic_m m_pos       = eq(min_cost,cost);
	    const unsigned long pos = bsf64(m_pos);	    

	    assert(pos < 15);

	    bestCost     = cost[pos];
	    bestSplit    = pos+1;
	    bestSplitDim = dim;	    
	    bestNumLeft  = lnum[pos];

	  }
      };


    if (unlikely(bestSplit == -1)) goto createLeaf;
    
    /* --------------------------- */
    /* --- sort left and right --- */
    /* --------------------------- */


    left.reset();
    right.reset();

    const unsigned int mid = partitionAABBs<L2_PREFETCH_DISTANCE>(aabb ,current.begin, current.end-1, bestSplit, bestSplitDim, centroidBoundsMin_2, scale, left, right);

  
    assert(current.begin + mid == current.begin + bestNumLeft);

    if (unlikely(current.begin + mid == current.begin || current.begin + mid == current.end)) 
      {
	std::cout << "warning: mid == current.begin || mid == current.end " << std::endl;
	DBG_PRINT(current);
	DBG_PRINT(mid);
	exit(0);
	goto createLeaf;
      }


    CHECK(
	  {
	    MIC_ALIGN Centroid_Scene_AABB test_left;
	    MIC_ALIGN Centroid_Scene_AABB test_right;
    
	    test_left.reset();
	    test_right.reset();

	    AABB * __restrict__ l   = aabb + current.begin;
	    AABB * __restrict__ r   = aabb + current.end;

	    unsigned int numLeft = 0;
	    unsigned int numRight = 0;
	    for (;l<r;l++)
	      {

		if (getBinSlot(l,bestSplitDim,centroidBoundsMin_2,scale) < bestSplit)
		  {
		    test_left.extend(*l);
		    numLeft++;

		  }
		else
		  {

		    test_right.extend(*l);	    
		    numRight++;
		  }
	      }

	    if (!(numLeft == bestNumLeft &&
		  test_left == left &&
		  test_right == right))
	      {
		PING;
		std::cout << "CHECKING SERIAL..." << std::endl;
		DBG_PRINT(current.items());
		DBG_PRINT(current.begin);
		DBG_PRINT(current.end);

		DBG_PRINT(bestNumLeft);
		DBG_PRINT(numLeft);
		DBG_PRINT(numRight);
		DBG_PRINT(test_left);
		DBG_PRINT(left);
		DBG_PRINT(test_right);
		DBG_PRINT(right);

		exit(0);
	      }


	    assert(numLeft == bestNumLeft);
	    assert(test_left == left);
	    assert(test_right == right);

	  }
	  );




    const unsigned int current_mid = current.begin + bestNumLeft;

    leftChild.init(left,current.begin,current_mid);
    rightChild.init(right,current_mid,current.end);

    //CHECK(PRINT(leftChild));
    assert(leftChild.sceneArea()     >= 0.0f);
    assert(leftChild.centroidArea()  >= 0.0f);
    assert(rightChild.sceneArea()    >= 0.0f);
    assert(rightChild.centroidArea() >= 0.0f);

    return true;
  }

  // =================================================================
  // =================================================================
  // =================================================================

  unsigned int ParallelQBVHBuilder::recurseSAH(BuildRecord &current, 
					       const unsigned int mode, 
					       const unsigned int threadID,
					       unsigned int &localNodeID,
					       unsigned int &localNodeIDs)
  {
    MIC_ALIGN BuildRecord record[4]; 
    MIC_ALIGN BuildRecord left, right;


    CHECK(
	  checkBuildRecord(current);
	  );

    const unsigned int coreID = threadID >> 2;

    unsigned int numSplits = 1;
    record[0] = current;

    left.prefetchL1Ex();
    right.prefetchL1Ex();

    while(1)
      {
	bool couldSplit = false;
	if (numSplits < 4)
	  {
	    int index = -1;
	    if (1) // mode == RECURSE || mode == PUSH_LOCAL)
	      {
		float maxArea = -1.0f;
		for (unsigned int i=0;i<numSplits;i++)
		  {
		    if (!record[i].isLeaf())
		      {
			assert(record[i].sceneArea() >= 0.0f);

			if (record[i].sceneArea() >= maxArea)
			  {
			    maxArea = record[i].sceneArea();
			    index = i;
			  }
		      }
		  }
	      }
	    else
	      {
		int max_items = 0;
		for (unsigned int i=0;i<numSplits;i++)
		  {
		    if (!record[i].isLeaf())
		      {
			if (record[i].items() > max_items)
			  {
			    max_items = record[i].items();
			    index = i;
			  }
		      }
		  }
	      } 

	    assert(index < 4);
	    if (index != -1)
	      {
		bool s = false;
		if (mode == PUSH_GLOBAL)
		  s = split_parallel(record[index],left,right);
		else
		  s = split(record[index],left,right);

		assert(numSplits > 0);
		record[numSplits].prefetchL1Ex();

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
		    continue;
		  }
	      }
	  }
	// ==================================================
	if (!couldSplit) break;
      }





    // could not generate any split, propagate leaf back to parent
    if ( numSplits == 1 )
      {
	if(current.items() <= 4)
	  {
	    current.createLeaf();

	    const unsigned int items  = current.items();
	    const unsigned int offset = current.begin;
	    assert(items  != 0);

	    globalNodePtr[current.nodeID].createLeaf(offset,items);
	    return numSplits;
	  }
	else 
	  {
	    numSplits = handleNonSplitableLeaves(current,record);
	  }      
      }

    assert(numSplits >= 2 && numSplits <= 4);

    // ==== aquire next four nodes ====
    const unsigned int currentIndex = getAtomicID(mode,localNodeID,localNodeIDs,4);
    if (unlikely(currentIndex >= max_nodes))
      {
	PRINT(currentIndex);
	PRINT(max_nodes);
	FATAL("not enough nodes allocated");
      }

    current.childrenID = currentIndex;

    // ==== init all four nodes ====
    const mic_f init_node = upconv16f((float*)ParallelQBVHBuilder::initQBVHNode);
    store16f_ngo(&globalNodePtr[currentIndex+0],init_node);
    store16f_ngo(&globalNodePtr[currentIndex+2],init_node);


    unsigned int nodes_to_process = 0;
    unsigned int node_to_process[4];

    for (unsigned int i=0;i<numSplits;i++)
      {
	record[i].cs_AABB.storeSceneAABB((float*)&globalNodePtr[currentIndex+i]);

	if (!record[i].isLeaf())
	  {
	    record[i].nodeID = currentIndex+i;
	    node_to_process[nodes_to_process++] = i;
	  }
	else
	  {
	    const unsigned int items  = record[i].items();
	    const unsigned int offset = record[i].begin;
	    globalNodePtr[currentIndex+i].createLeaf(offset,items);
	  }
      
      }
    globalNodePtr[current.nodeID].createNode(current.childrenID,numSplits);

    assert(current.childrenID > 0);
    assert(numSplits > 0);


    if (nodes_to_process > 0)
      {
	for (unsigned int i=0;i<nodes_to_process;i++)
	  {
	    const unsigned int index = node_to_process[i];
	    if (mode == PUSH_GLOBAL)
	      {
		atomicBuildRecordStack.push_nolock(record[index]);		
	      }
	    else if (mode == PUSH_LOCAL)
	      {
		if (localAtomicBuildRecordStack[coreID].usedSlots() < NUM_LOCAL_ATOMIC_BUILD_RECORDS-1 && 
		    record[index].items() > LOCAL_BUILD_RECORD_STACK_PUSH_THRESHOLD )
		  {
		    while (localAtomicBuildRecordStack[coreID].push(record[index]) == false)
		      {
			delayThread(1024);
			CHECK(std::cout << "LOCAL ENQUEUE ERROR " << coreID << std::endl);
		      }
		  }
		else
		  {
		    recurseSAH(record[node_to_process[i]],mode,threadID,localNodeID,localNodeIDs);
		  }

	      }
	    else
	      {
		recurseSAH(record[node_to_process[i]],mode,threadID,localNodeID,localNodeIDs);
	      }
	  }
      }

    return numSplits;  
  }



  // =================================================================
  // =================================================================
  // =================================================================

  void ParallelQBVHBuilder::thread_convertToOptimizedQBVH(const unsigned int threadID)
  {
    if (unlikely(threadID >= splitData.worker_threads)) return;

    const unsigned int WORKERS  = splitData.worker_threads;
    const unsigned int workerID = threadID;
    const unsigned int nodes =  qbvh_nodes;
    const unsigned int primitives_per_block = nodes / WORKERS;
    const unsigned int startID = threadID * primitives_per_block;
    const unsigned int endID  = (workerID != WORKERS-1) ? startID + primitives_per_block : nodes;

    QBVHNode *qbvh = (QBVHNode*)globalNodePtr;

    QBVHNode * __restrict__  qptr = qbvh          + startID;
    BVHNode  * __restrict__  bptr = globalNodePtr + startID*4;

    for (unsigned int n=startID;n<endID;n++,qptr++,bptr+=4)
      {
	prefetch<PFHINT_L1EX>(qptr  );
	prefetch<PFHINT_L1EX>(qptr+1);
	prefetch<PFHINT_L1>(bptr+4);
	prefetch<PFHINT_L2EX>(qptr+8);
	prefetch<PFHINT_L2>(bptr+8);

	QBVHNode tmp;
	for (int i=0;i<4;i++)
	  {
	    tmp.m_min[i].x = bptr[i].m_min[0];
	    tmp.m_min[i].y = bptr[i].m_min[1];
	    tmp.m_min[i].z = bptr[i].m_min[2];

	    tmp.m_max[i].x = bptr[i].m_max[0];
	    tmp.m_max[i].y = bptr[i].m_max[1];
	    tmp.m_max[i].z = bptr[i].m_max[2];
	    tmp.m_max[i].d = bptr[i].ext_max.t;

	    if (!bvhLeaf(bptr[i].ext_min.t))
	      {
		tmp.m_min[i].d = qbvhCreateNode(bvhChildID(bptr[i].ext_min.t)>>2,0); // bvhChildren(bptr[i].ext_min.t)
	      }
	    else
	      {
		tmp.m_min[i].d = (bptr[i].ext_min.t ^ BVH_LEAF_MASK) | QBVH_LEAF_MASK;
	      }	  
	  }
	*qptr = tmp;
      }
  
  }


  void ParallelQBVHBuilder::convertToOptimizedQBVH()
  {
    QBVHTaskScheduler::dispatchTask(ParallelQBVHBuilder::thread_convertToOptimizedQBVH);
  }

  // =================================================================
  // =================================================================
  // =================================================================

  void ParallelQBVHBuilder::thread_createTriangleAccel(const unsigned int threadID)
  {
    if (unlikely(threadID >= splitData.worker_threads)) return;

    const unsigned int WORKERS  = splitData.worker_threads;
    const unsigned int workerID = threadID;
    const unsigned int primitives_per_block = triangles / WORKERS;
    const unsigned int startID = threadID * primitives_per_block;
    const unsigned int endID  = (workerID != WORKERS-1) ? startID + primitives_per_block : triangles;

    TriangleAccel * __restrict__  acc = globalAccelPtr + startID;
    const AABB    * __restrict__  bptr = aabb + startID;

    const Vec3fa          * __restrict__ const vptr = globalVertexPtr;
    const BuilderTriangle * __restrict__       tptr = globalBuildTrianglePtr;


    prefetch<PFHINT_L1>(bptr);
    for (unsigned int j=startID;j<endID;j++,bptr++,acc++)
      {
	prefetch<PFHINT_L1>(bptr + 2);
	prefetch<PFHINT_L2>(bptr + L2_PREFETCH_DISTANCE);

	const int index    = bptr[0].ext_min.t;
	const int index_pf = bptr[1].ext_min.t;

	prefetch<PFHINT_NT>(&tptr[index_pf]);
	prefetch<PFHINT_NT>(&tptr[index]);

	const mic_f tAccel = computeTriangleAccel(index,tptr,vptr);

	store16f_ngo(acc,tAccel);
      }

  }


  void ParallelQBVHBuilder::computeSAHfromReduceBins()
  {
    Bin16 &bin16 = threadBin16[0];

    const unsigned int items = splitData.rec.items();
    const float voxelArea = splitData.rec.cs_AABB.sceneArea();
    const mic_f centroidDiagonal_2  = splitData.rec.cs_AABB.centroidDiagonal() * 2.0f;
    const mic_f centroidBoundsMin_2 = splitData.rec.cs_AABB.minCentroidAABB()  * 2.0f;
    const mic_f scale = sel(nz(centroidDiagonal_2),rcp(centroidDiagonal_2) * mic_f(16.0f * 0.99f),mic_f::zero());

    BestSplitData bestSplit;
    bestSplit.init(items * voxelArea);

    for (unsigned int dim = 0;dim < 3;dim++) 
      {
	if (unlikely(centroidDiagonal_2[dim] == 0.0f)) continue;
      
	const mic_f rArea   = bin16.prefix_area_rl(dim);
	const mic_f lArea   = bin16.prefix_area_lr(dim);      
	const mic_i lnum    = bin16.prefix_count(dim);
      
	const mic_i rnum    = mic_i(items) - lnum;
	const mic_i lblocks = (lnum + mic_i(3)) >> 2;
	const mic_i rblocks = (rnum + mic_i(3)) >> 2;
	const mic_m m_lnum  = eqz(lnum);
	const mic_m m_rnum  = eqz(rnum);
	const mic_f cost    = sel(m_lnum|m_rnum,mic_f::inf(),lArea * mic_f(lblocks) + rArea * mic_f(rblocks) + voxelArea );
	if (lt(cost,mic_f(bestSplit.cost)))
	  {
	    const mic_f min_cost    = set_min16(cost); 
	    const mic_m m_pos       = eq(min_cost,cost);
	    const unsigned long pos = bsf64(m_pos);	    
	      
	    assert(pos < 15);

	    bestSplit.cost  = cost[pos];
	    bestSplit.index = pos+1;
	    bestSplit.dim   = dim;	    
	    bestSplit.left  = lnum[pos];
	  }
      };


    splitData.rec.splitIndex       = bestSplit.index;
    splitData.rec.splitNumLeft     = bestSplit.left;
    splitData.bestSplit            = bestSplit.index;
    splitData.bestSplitDim         = bestSplit.dim;
    splitData.centroidBoundsMin_2  = centroidBoundsMin_2;
    splitData.scale                = scale;

  }

  // ==============================================================================================================
  // ==============================================================================================================
  // ==============================================================================================================


  void ParallelQBVHBuilder::fastbinning(const unsigned int start, const unsigned int end, Bin16 &bin16, const unsigned int threadID)
  {
    unsigned int thread_start = start;
    unsigned int thread_end   = end;

    prefetch<PFHINT_NT>(&aabb[thread_start]);
    prefetch<PFHINT_NT>(&aabb[thread_end-1]);

    prefetch<PFHINT_L1EX>(&tmp_aabb[thread_start]);
    prefetch<PFHINT_L1EX>(&tmp_aabb[thread_end-1]);

    const mic_f init_min = mic_f::inf();
    const mic_f init_max = mic_f::minus_inf();
    const mic_i zero     = mic_i::zero();

    mic_f min_x0,min_x1,min_x2;
    mic_f min_y0,min_y1,min_y2;
    mic_f min_z0,min_z1,min_z2;
    mic_f max_x0,max_x1,max_x2;
    mic_f max_y0,max_y1,max_y2;
    mic_f max_z0,max_z1,max_z2;
    mic_i count0,count1,count2;

    min_x0 = init_min;
    min_x1 = init_min;
    min_x2 = init_min;
    min_y0 = init_min;
    min_y1 = init_min;
    min_y2 = init_min;
    min_z0 = init_min;
    min_z1 = init_min;
    min_z2 = init_min;

    max_x0 = init_max;
    max_x1 = init_max;
    max_x2 = init_max;
    max_y0 = init_max;
    max_y1 = init_max;
    max_y2 = init_max;
    max_z0 = init_max;
    max_z1 = init_max;
    max_z2 = init_max;

    count0 = zero;
    count1 = zero;
    count2 = zero;

    const mic_f centroidBoundsMin_2 = splitData.rec.cs_AABB.minCentroidAABB()  * 2.0f;
    const mic_f centroidDiagonal_2  = splitData.rec.cs_AABB.centroidDiagonal() * 2.0f;
    const mic_f scale = sel(nz(centroidDiagonal_2),rcp(centroidDiagonal_2) * mic_f(16.0f * 0.99f),mic_f::zero());

    const mic_i id = mic_i::identity();


    if (likely(thread_start % 2 != 0))
      {
	const mic_f b_min = upconv4f((float*)&aabb[thread_start].m_min);
	const mic_f b_max = upconv4f((float*)&aabb[thread_start].m_max);    

	const mic_f centroid_2 = b_min + b_max;
	const mic_i binID = mic_i((centroid_2 - centroidBoundsMin_2)*scale);

	assert(0 <= binID[0] && binID[0] < 16);
	assert(0 <= binID[1] && binID[1] < 16);
	assert(0 <= binID[2] && binID[2] < 16);

	const mic_m m_update_x = eq(id,swAAAA_i(binID));
	const mic_m m_update_y = eq(id,swBBBB_i(binID));
	const mic_m m_update_z = eq(id,swCCCC_i(binID));

	min_x0 = mask_min(m_update_x,min_x0,min_x0,swAAAA(b_min));
	min_y0 = mask_min(m_update_x,min_y0,min_y0,swBBBB(b_min));
	min_z0 = mask_min(m_update_x,min_z0,min_z0,swCCCC(b_min));      
	// ------------------------------------------------------------------------
	max_x0 = mask_max(m_update_x,max_x0,max_x0,swAAAA(b_max));
	max_y0 = mask_max(m_update_x,max_y0,max_y0,swBBBB(b_max));
	max_z0 = mask_max(m_update_x,max_z0,max_z0,swCCCC(b_max));
	// ------------------------------------------------------------------------
	min_x1 = mask_min(m_update_y,min_x1,min_x1,swAAAA(b_min));
	min_y1 = mask_min(m_update_y,min_y1,min_y1,swBBBB(b_min));
	min_z1 = mask_min(m_update_y,min_z1,min_z1,swCCCC(b_min));
	// ------------------------------------------------------------------------      
	max_x1 = mask_max(m_update_y,max_x1,max_x1,swAAAA(b_max));
	max_y1 = mask_max(m_update_y,max_y1,max_y1,swBBBB(b_max));
	max_z1 = mask_max(m_update_y,max_z1,max_z1,swCCCC(b_max));
	// ------------------------------------------------------------------------
	min_x2 = mask_min(m_update_z,min_x2,min_x2,swAAAA(b_min));
	min_y2 = mask_min(m_update_z,min_y2,min_y2,swBBBB(b_min));
	min_z2 = mask_min(m_update_z,min_z2,min_z2,swCCCC(b_min));
	// ------------------------------------------------------------------------      
	max_x2 = mask_max(m_update_z,max_x2,max_x2,swAAAA(b_max));
	max_y2 = mask_max(m_update_z,max_y2,max_y2,swBBBB(b_max));
	max_z2 = mask_max(m_update_z,max_z2,max_z2,swCCCC(b_max));
	// ------------------------------------------------------------------------
	count0 = mask_add(m_update_x,count0,count0,mic_i::one());
	count1 = mask_add(m_update_y,count1,count1,mic_i::one());
	count2 = mask_add(m_update_z,count2,count2,mic_i::one());      
	// ------------------------------------------------------------------------
	tmp_aabb[thread_start] = aabb[thread_start];
	thread_start++;
      }

    if (likely(thread_end % 2 != 0))
      {
	const mic_f b_min = upconv4f((float*)&aabb[thread_end-1].m_min);
	const mic_f b_max = upconv4f((float*)&aabb[thread_end-1].m_max);    

	const mic_f centroid_2 = b_min + b_max;
	const mic_i binID = mic_i((centroid_2 - centroidBoundsMin_2)*scale);

	assert(0 <= binID[0] && binID[0] < 16);
	assert(0 <= binID[1] && binID[1] < 16);
	assert(0 <= binID[2] && binID[2] < 16);

	const mic_m m_update_x = eq(id,swAAAA_i(binID));
	const mic_m m_update_y = eq(id,swBBBB_i(binID));
	const mic_m m_update_z = eq(id,swCCCC_i(binID));

	min_x0 = mask_min(m_update_x,min_x0,min_x0,swAAAA(b_min));
	min_y0 = mask_min(m_update_x,min_y0,min_y0,swBBBB(b_min));
	min_z0 = mask_min(m_update_x,min_z0,min_z0,swCCCC(b_min));      
	// ------------------------------------------------------------------------
	max_x0 = mask_max(m_update_x,max_x0,max_x0,swAAAA(b_max));
	max_y0 = mask_max(m_update_x,max_y0,max_y0,swBBBB(b_max));
	max_z0 = mask_max(m_update_x,max_z0,max_z0,swCCCC(b_max));
	// ------------------------------------------------------------------------
	min_x1 = mask_min(m_update_y,min_x1,min_x1,swAAAA(b_min));
	min_y1 = mask_min(m_update_y,min_y1,min_y1,swBBBB(b_min));
	min_z1 = mask_min(m_update_y,min_z1,min_z1,swCCCC(b_min));
	// ------------------------------------------------------------------------      
	max_x1 = mask_max(m_update_y,max_x1,max_x1,swAAAA(b_max));
	max_y1 = mask_max(m_update_y,max_y1,max_y1,swBBBB(b_max));
	max_z1 = mask_max(m_update_y,max_z1,max_z1,swCCCC(b_max));
	// ------------------------------------------------------------------------
	min_x2 = mask_min(m_update_z,min_x2,min_x2,swAAAA(b_min));
	min_y2 = mask_min(m_update_z,min_y2,min_y2,swBBBB(b_min));
	min_z2 = mask_min(m_update_z,min_z2,min_z2,swCCCC(b_min));
	// ------------------------------------------------------------------------      
	max_x2 = mask_max(m_update_z,max_x2,max_x2,swAAAA(b_max));
	max_y2 = mask_max(m_update_z,max_y2,max_y2,swBBBB(b_max));
	max_z2 = mask_max(m_update_z,max_z2,max_z2,swCCCC(b_max));
	// ------------------------------------------------------------------------
	count0 = mask_add(m_update_x,count0,count0,mic_i::one());
	count1 = mask_add(m_update_y,count1,count1,mic_i::one());
	count2 = mask_add(m_update_z,count2,count2,mic_i::one());      
	// ------------------------------------------------------------------------
	tmp_aabb[thread_end-1] = aabb[thread_end-1];
	thread_end--;
      }



    const AABB * __restrict__ aptr = aabb + thread_start;
    AABB * __restrict__ tmp_aptr   = tmp_aabb + thread_start;

    prefetch<PFHINT_NT>(aptr);
    prefetch<PFHINT_L2>(aptr+2);
    prefetch<PFHINT_L2>(aptr+4);
    prefetch<PFHINT_L2>(aptr+6);
    prefetch<PFHINT_L2>(aptr+8);

    for (unsigned int j = thread_start;j < thread_end;j+=2,aptr+=2,tmp_aptr+=2)
      {
	prefetch<PFHINT_NT>(aptr+2);
	prefetch<PFHINT_L2>(aptr+8);

	assert((unsigned long)aptr % 64 == 0);
	assert((unsigned long)tmp_aptr % 64 == 0);

	store16f_ngo(tmp_aptr,*(mic_f*)aptr);

#pragma unroll
	for (int i=0;i<2;i++)
	  {
	    const mic_f b_min = upconv4f((float*)&aptr[i].m_min);
	    const mic_f b_max = upconv4f((float*)&aptr[i].m_max);    

	    const mic_f centroid_2 = b_min + b_max;
	    const mic_i binID = mic_i((centroid_2 - centroidBoundsMin_2)*scale);

	    assert(0 <= binID[0] && binID[0] < 16);
	    assert(0 <= binID[1] && binID[1] < 16);
	    assert(0 <= binID[2] && binID[2] < 16);

	    const mic_m m_update_x = eq(id,swAAAA_i(binID));
	    const mic_m m_update_y = eq(id,swBBBB_i(binID));
	    const mic_m m_update_z = eq(id,swCCCC_i(binID));

	    min_x0 = mask_min(m_update_x,min_x0,min_x0,swAAAA(b_min));
	    min_y0 = mask_min(m_update_x,min_y0,min_y0,swBBBB(b_min));
	    min_z0 = mask_min(m_update_x,min_z0,min_z0,swCCCC(b_min));
      
	    max_x0 = mask_max(m_update_x,max_x0,max_x0,swAAAA(b_max));
	    max_y0 = mask_max(m_update_x,max_y0,max_y0,swBBBB(b_max));
	    max_z0 = mask_max(m_update_x,max_z0,max_z0,swCCCC(b_max));
	    // ------------------------------------------------------------------------
	    min_x1 = mask_min(m_update_y,min_x1,min_x1,swAAAA(b_min));
	    min_y1 = mask_min(m_update_y,min_y1,min_y1,swBBBB(b_min));
	    min_z1 = mask_min(m_update_y,min_z1,min_z1,swCCCC(b_min));
      
	    max_x1 = mask_max(m_update_y,max_x1,max_x1,swAAAA(b_max));
	    max_y1 = mask_max(m_update_y,max_y1,max_y1,swBBBB(b_max));
	    max_z1 = mask_max(m_update_y,max_z1,max_z1,swCCCC(b_max));
	    // ------------------------------------------------------------------------
	    min_x2 = mask_min(m_update_z,min_x2,min_x2,swAAAA(b_min));
	    min_y2 = mask_min(m_update_z,min_y2,min_y2,swBBBB(b_min));
	    min_z2 = mask_min(m_update_z,min_z2,min_z2,swCCCC(b_min));
      
	    max_x2 = mask_max(m_update_z,max_x2,max_x2,swAAAA(b_max));
	    max_y2 = mask_max(m_update_z,max_y2,max_y2,swBBBB(b_max));
	    max_z2 = mask_max(m_update_z,max_z2,max_z2,swCCCC(b_max));
      
	    count0 = mask_add(m_update_x,count0,count0,mic_i::one());
	    count1 = mask_add(m_update_y,count1,count1,mic_i::one());
	    count2 = mask_add(m_update_z,count2,count2,mic_i::one());      

	  }
      }

    bin16.prefetchL1EX();

    bin16.min_x[0] = min_x0;
    bin16.min_y[0] = min_y0;
    bin16.min_z[0] = min_z0;
    bin16.max_x[0] = max_x0;
    bin16.max_y[0] = max_y0;
    bin16.max_z[0] = max_z0;

    bin16.min_x[1] = min_x1;
    bin16.min_y[1] = min_y1;
    bin16.min_z[1] = min_z1;
    bin16.max_x[1] = max_x1;
    bin16.max_y[1] = max_y1;
    bin16.max_z[1] = max_z1;

    bin16.min_x[2] = min_x2;
    bin16.min_y[2] = min_y2;
    bin16.min_z[2] = min_z2;
    bin16.max_x[2] = max_x2;
    bin16.max_y[2] = max_y2;
    bin16.max_z[2] = max_z2;

    bin16.count[0] = count0;
    bin16.count[1] = count1;
    bin16.count[2] = count2;

    bin16.thread_count[0] = count0; 
    bin16.thread_count[1] = count1; 
    bin16.thread_count[2] = count2; 

  }


  // ==================================================================================================================
  // ==================================================================================================================
  // ==================================================================================================================


  void ParallelQBVHBuilder::thread_partition_parallel_tmp_copy(const unsigned int threadID)
  {
    splitData.prefetchL2();
    if (unlikely(threadID >= splitData.worker_threads)) return;

    const unsigned int WORKERS      = splitData.worker_threads;
    const unsigned int workerID     = threadID;
    const unsigned int items        = splitData.rec.items();

    const unsigned int bestSplit     = splitData.bestSplit;
    const unsigned int bestSplitDim  = splitData.bestSplitDim;
    const mic_f centroidBoundsMin_2  = splitData.centroidBoundsMin_2;
    const mic_f scale                = splitData.scale;

    const unsigned int primitives_per_thread = splitData.work_items_per_thread;
    const unsigned int thread_begin = splitData.rec.begin + threadID * primitives_per_thread;
    const unsigned int thread_end   = (workerID != WORKERS-1) ? thread_begin + primitives_per_thread : splitData.rec.end;

    const mic_f c = mic_f(centroidBoundsMin_2[bestSplitDim]);
    const mic_f s = mic_f(scale[bestSplitDim]);


    // compute items per thread that go to the 'left' and to the 'right'
    const mic_i lnum    = prefix_sum(threadBin16[threadID].thread_count[splitData.bestSplitDim]);
    const unsigned int local_numLeft = lnum[splitData.bestSplit-1];
    const unsigned int local_numRight = (thread_end-thread_begin) - lnum[splitData.bestSplit-1];
 
    const unsigned int thread_start_left  = lCounter.add(local_numLeft);
    const unsigned int thread_start_right = rCounter.add(local_numRight);



    MIC_ALIGN Centroid_Scene_AABB local_left;
    MIC_ALIGN Centroid_Scene_AABB local_right;
    
    local_left.reset();
    local_right.reset();

    AABB * __restrict__ l_source = tmp_aabb + thread_begin;
    AABB * __restrict__ r_source = tmp_aabb + thread_end;

    AABB * __restrict__ l_dest     = aabb + splitData.rec.begin + thread_start_left;
    AABB * __restrict__ r_dest     = aabb + splitData.rec.begin + thread_start_right + splitData.rec.splitNumLeft;

    for (;l_source<r_source;)
      {
	prefetch<PFHINT_NT>(l_source+2);
	prefetch<PFHINT_L2>(l_source + L2_PREFETCH_DISTANCE);

	if (likely(lt_split(l_source,bestSplitDim,c,s,mic_f(bestSplit))))
	  {
	    prefetch<PFHINT_L1EX>(l_dest+2);
	    prefetch<PFHINT_L2EX>(l_dest + L2_PREFETCH_DISTANCE);
	    local_left.extend(*l_source); 
	    *l_dest++ = *l_source++;
	    evictL1(l_dest-2);
	    evictL1(l_source-2);
	  }
	else
	  {
	    prefetch<PFHINT_L1EX>(r_dest+2);
	    prefetch<PFHINT_L2EX>(r_dest + L2_PREFETCH_DISTANCE);
	    local_right.extend(*l_source); 
	    *r_dest++ = *l_source++;
	    evictL1(r_dest-2);
	    evictL1(l_source-2);

	  }
      }

    splitData.left.extend_atomic(local_left); 
    splitData.right.extend_atomic(local_right);  
  }


  // =============================================================================================
  // =============================================================================================
  // =============================================================================================



  void ParallelQBVHBuilder::reduceBin16(const unsigned int currentThreadID,
					const unsigned int childThreadID)
  {
    if (unlikely(childThreadID >= splitData.worker_threads)) return;

    Bin16  * __restrict__ b0 = ParallelQBVHBuilder::getThreadBin16Ptr(currentThreadID);
    Bin16  * __restrict__ b1 = ParallelQBVHBuilder::getThreadBin16Ptr(childThreadID);
    b0->merge_with_prefetch(*b1);
  }

  void ParallelQBVHBuilder::thread_binning(const unsigned int threadID)
  {
    splitData.prefetchL2();

    const unsigned int WORKERS  = splitData.worker_threads;
    const unsigned int workerID = threadID;       
    Bin16 &bin16 = ParallelQBVHBuilder::threadBin16[threadID];

    if (likely(threadID < splitData.worker_threads)) 
      {

	const float voxelArea              = splitData.rec.cs_AABB.sceneArea();
	const unsigned int items           = splitData.rec.items();
	const unsigned int items_per_block = splitData.work_items_per_thread; 
	const unsigned int startID         = splitData.rec.begin + workerID * items_per_block;
	const unsigned int endID           = (workerID != WORKERS-1) ? startID + items_per_block : splitData.rec.end;
  
	assert(endID - startID > 0);


	const mic_f centroidDiagonal_2     = splitData.rec.cs_AABB.centroidDiagonal() * 2.0f;
	const mic_f centroidBoundsMin_2    = splitData.rec.cs_AABB.minCentroidAABB()  * 2.0f;
	const mic_f scale                  = sel(nz(centroidDiagonal_2),rcp(centroidDiagonal_2) * mic_f(16.0f * 0.99f),mic_f::zero());

	fastbinning(startID,endID,bin16,threadID);
      }
    taskBarrier.syncWithReduction(workerID,NUM_TOTAL_THREADS,&reduceBin16);
  }


  bool ParallelQBVHBuilder::split_parallel(BuildRecord &current,
					   BuildRecord &leftChild,
					   BuildRecord &rightChild)
  {
    CHECK( PING );
    const unsigned int items = current.end - current.begin;
    const float voxelArea    = current.cs_AABB.sceneArea();
    const float centroidArea = current.cs_AABB.centroidArea();
    assert(items > 0);
  
    if (items <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD  || 
	voxelArea <= 0.0f ||
	eqz(0x7777,current.cs_AABB.centroidDiagonal()) == 0x7777
	) 
      {
	current.createLeaf();
	return false;
      }

    /* --------------------- */
    /* --- bin centroids --- */
    /* --------------------- */


#if 0
    Timer timer; 
    timer.start();
#endif
    bool s = true;


    if (unlikely(current.items() <= BUILD_RECORD_PARALLEL_SPLIT_THRESHOLD))
      {
	s = split(current,leftChild,rightChild);	
      }
    else
      {
	splitData.worker_threads = 4 * min( NUM_TOTAL_CORES, current.items() / BUILD_RECORD_PARALLEL_SPLIT_THRESHOLD );

	if (splitData.worker_threads % 4) splitData.worker_threads -= splitData.worker_threads % 4;

	if (unlikely(splitData.worker_threads < 4)) FATAL("HERE");

	splitData.work_items_per_thread = current.items() / splitData.worker_threads;
	splitData.rec = current;
	splitData.rec.splitIndex   = -1;
	splitData.rec.splitNumLeft = -1;
	splitData.left.reset();
	splitData.right.reset();


	QBVHTaskScheduler::dispatchTask( ParallelQBVHBuilder::thread_binning );


	ParallelQBVHBuilder::computeSAHfromReduceBins();

	if (unlikely(splitData.rec.splitIndex == -1)) 
	  {
	    splitData.worker_threads = NUM_TOTAL_THREADS;
	    current.createLeaf();
	    return false;
	  }



	const unsigned int mid = partition_parallel(splitData.rec,
						    splitData.bestSplit,
						    splitData.bestSplitDim,
						    splitData.centroidBoundsMin_2,
						    splitData.scale,
						    splitData.left,
						    splitData.right);



	assert(mid == splitData.rec.begin + splitData.rec.splitNumLeft);

	//const unsigned int mid = current.begin + splitData.rec.splitNumLeft;

	CHECK(
	      AABB left;
	      splitData.left.storeSceneAABB((float*)&left);
	      DBG_PRINT(left);
	      for (unsigned int i=current.begin;i<mid;i++) 
		{
		  if (!left.enclose(aabb[i]))
		    {
		      DBG_PRINT(i);
		      DBG_PRINT(aabb[i]);
		    }

		  assert(left.enclose(aabb[i]));        
		}
	      );
	CHECK(
	      AABB right;
	      splitData.right.storeSceneAABB((float*)&right);
	      DBG_PRINT(right);

	      for (unsigned int i=mid;i<current.end;i++) 
		assert(right.enclose(aabb[i])); 
	      );

	assert(mid != (unsigned int)-1);

	leftChild.init(splitData.left,current.begin,mid);
	rightChild.init(splitData.right,mid,current.end);

	assert(leftChild.sceneArea()     > 0.0f);
	assert(leftChild.centroidArea()  >= 0.0f);
	assert(rightChild.sceneArea()    > 0.0f);
	assert(rightChild.centroidArea() >= 0.0f);

	CHECK(
	      checkBuildRecord(leftChild);
	      checkBuildRecord(rightChild);
	      );

      }


    splitData.worker_threads = NUM_TOTAL_THREADS;

    return s;
  }

  // =============================================================================================
  // =============================================================================================
  // =============================================================================================

  unsigned int ParallelQBVHBuilder::partition_parallel(const BuildRecord &current,
						       const unsigned int bestSplit,
						       const unsigned int bestSplitDim,
						       const mic_f &centroidBoundsMin_2,
						       const mic_f &scale,
						       Centroid_Scene_AABB &left,
						       Centroid_Scene_AABB &right)
  {
    splitData.rec = current;
    splitData.left.reset();
    splitData.right.reset();
    splitData.bestSplit = bestSplit;
    splitData.bestSplitDim = bestSplitDim;
    splitData.centroidBoundsMin_2 = centroidBoundsMin_2;
    splitData.scale = scale;


    lCounter.prefetchEx();
    rCounter.prefetchEx();
    lCounter.reset(0);
    rCounter.reset(0); 

    dispatchTask(ParallelQBVHBuilder::thread_partition_parallel_tmp_copy);
    assert(lCounter.val() == splitData.rec.splitNumLeft);
    assert(rCounter.val() == current.items() -  splitData.rec.splitNumLeft);


    const unsigned int splitIndex = current.begin + splitData.rec.splitNumLeft;

    CHECK(
	  MIC_ALIGN Centroid_Scene_AABB test_left;
	  MIC_ALIGN Centroid_Scene_AABB test_right;
    
	  test_left.reset();
	  test_right.reset();

	  AABB * __restrict__ l   = aabb + current.begin;
	  AABB * __restrict__ r   = aabb + current.end;

	  DBG_PRINT(current.items());
	  DBG_PRINT(current.begin);
	  DBG_PRINT(current.end);
	  DBG_PRINT(splitIndex);

	  unsigned int numLeft = 0;
	  unsigned int numRight = 0;
	  for (;l<r;l++)
	    {

	      if (getBinSlot(l,bestSplitDim,centroidBoundsMin_2,scale) < bestSplit)
		{
		  test_left.extend(*l);
		  numLeft++;

		}
	      else
		{

		  test_right.extend(*l);	    
		  numRight++;
		}
	    }

	  if (!(numLeft == splitData.rec.splitNumLeft &&
		test_left == splitData.left &&
		test_right == splitData.right))
	    {
	      std::cout << "CHECKING PARTITION PARALLEL..." << std::endl;
	      DBG_PRINT(splitData.rec.splitNumLeft);
		
	      DBG_PRINT(numLeft);
	      DBG_PRINT(numRight);
	      DBG_PRINT(test_left);
	      DBG_PRINT(test_right);
	      DBG_PRINT(splitData.left);
	      DBG_PRINT(splitData.right);
	      exit(0);
	    }

	  assert(numLeft == splitData.rec.splitNumLeft);
	  assert(test_left == splitData.left);
	  assert(test_right == splitData.right);

	  );

    CHECK(
	  const bool checkLeft = checkRange(current.begin,splitIndex,0);	
	  if (checkLeft != true) 
	    { FATAL("CHECKLEFT"); }

	  const bool checkRight = checkRange(splitIndex,current.end,1); 
	  if (checkRight != true) 
	    { FATAL("CHECKRIGHT"); }
	  );

    return splitIndex;

  }


  unsigned int ParallelQBVHBuilder::handleNonSplitableLeaves(BuildRecord &current, 
							     BuildRecord record[4])
  {
    unsigned int numSplits = 1;
    CHECK(std::cout << "too many items in leaf -> brute force splitting: " << current.items() << std::endl);
    const unsigned int items  = current.items();

    for (int i=0;i<4;i++)
      {
	record[i] = current;
	record[i].createLeaf();
      }

    record[0].end   = record[0].begin+4;
    record[1].begin = record[0].end;


    if (current.items() <= 8)
      {
	record[1].end   = record[1].begin + items - 4;
	numSplits = 2;
      }
    else if (current.items() <= 12)
      {
	record[1].end   = record[1].begin+4;
	record[2].begin = record[1].end;
	record[2].end   = record[2].begin + items - 8;
	numSplits = 3;
      }
    else if (current.items() <= 16)
      {
	record[1].end   = record[1].begin+4;
	record[2].begin = record[1].end;
	record[2].end   = record[2].begin+4;
	record[3].begin = record[2].end;
	record[3].end   = record[3].begin + items - 12;
	numSplits = 4;
      }
    else // current.items() > 16, should be very rare
      {
	record[1].end   = record[1].begin+4;
	record[2].begin = record[1].end;
	record[2].end   = record[2].begin+4;
	numSplits = 4;
	record[3].createNode();
	record[3].begin = record[2].begin+4;
	//record[3].end   = record[2].begin+8;
	
	CHECK(std::cout << "too many items in leaf:" << current.items() << std::endl);
      }


    return numSplits;  
  }


  // =================================================================
  // =================================================================
  // =================================================================

  void ParallelQBVHBuilder::transferBuildRecordsFromGlobalToLocalStack()
  {
    CHECK( 
	  
	  DBG_PRINT(numActiveLocalWorkQueues) 
	   );
    for (unsigned int i=0;i<NUM_TOTAL_CORES;i++)
      {
	atomicBuildRecordStack.get_addr(i*2+0)->prefetchL2Ex();
	atomicBuildRecordStack.get_addr(i*2+1)->prefetchL2Ex();

	mic_f *ptr = (mic_f*)localAtomicBuildRecordStack[i].get_addr(0);
	prefetch<PFHINT_L2EX>(ptr+0);
	prefetch<PFHINT_L2EX>(ptr+1);
	localAtomicBuildRecordStack[i].init();
      }

    unsigned int slots = min(numActiveLocalWorkQueues,atomicBuildRecordStack.usedSlots());

    unsigned int distribution[64];
    for (unsigned int i=0;i<slots;i++)
      {
	distribution[i] = atomicBuildRecordStack.t[i].items();
	CHECK( checkBuildRecord( atomicBuildRecordStack.get(i) ) );
	localAtomicBuildRecordStack[i].push_nolock(atomicBuildRecordStack.get(i));
      }
    for (unsigned int i=slots;i<atomicBuildRecordStack.usedSlots();i++)
      {
	unsigned int smallest = distribution[0];
	unsigned int smallest_index = 0;
	for (unsigned int j=1;j<numActiveLocalWorkQueues;j++)
	  if (distribution[j] < smallest)
	    {
	      smallest = distribution[j];
	      smallest_index = j;
	    }
	distribution[smallest_index] += atomicBuildRecordStack.t[i].items();
	CHECK( checkBuildRecord( atomicBuildRecordStack.get(i) ) );
	localAtomicBuildRecordStack[smallest_index].push_nolock(atomicBuildRecordStack.get(i));
      }
  }



  void ParallelQBVHBuilder::build()
  {
    CHECK( PING );

    MIC_ALIGN BuildRecord record;
    TIMER(Timer timer; unsigned int cycles = 0);

    if (unlikely(!triangles)) {
      // should actually never happen
      const mic_f init_node = upconv16f((float*)ParallelQBVHBuilder::initQBVHNode);
      store16f_ngo(&globalNodePtr[0],init_node);
      store16f_ngo(&globalNodePtr[2],init_node);
      store16f_ngo(&globalNodePtr[4],init_node);
      store16f_ngo(&globalNodePtr[6],init_node);
      return;
    }
    
    assert(triangles >= 1);

    splitData.worker_threads = NUM_TOTAL_THREADS;

    CHECK(
	  DBG_PRINT(globalNodePtr);
	  DBG_PRINT(globalVertexPtr);
	  DBG_PRINT(globalBuildTrianglePtr);
	  DBG_PRINT(globalAccelPtr);
	  DBG_PRINT(numActiveLocalWorkQueues);
	  DBG_PRINT(NUM_TOTAL_THREADS);
	  DBG_PRINT(NUM_TOTAL_CORES);
	  DBG_PRINT(triangles);
	  DBG_PRINT(aabb);
	  );

    
    // ===============================================================
    TIMER(timer.start());
    global_cs_AABB.reset();
    QBVHTaskScheduler::dispatchTask( ParallelQBVHBuilder::thread_createAABBs );
    global_cs_AABB.storeSceneAABB((float*)&globalNodePtr[0]);
    TIMER(cycles = (unsigned long)timer.stop());    
    TIMER(std::cout << "compute AABBs in " << cycles << " cycles" << std::endl << std::flush);
    // ===============================================================

    CHECK(
	  DBG_PRINT(global_cs_AABB);
	  );


    recursePtr = recurseSAH;
    //while(1) 
      {
    ParallelQBVHBuilder::atomicBuildRecordStack.init();
    ParallelQBVHBuilder::atomicID.reset(4);

    record.init(global_cs_AABB, 0, triangles);
    record.nodeID = 0;

    atomicBuildRecordStack.push_nolock(record);  

    // ===============================================================
    TIMER(timer.start());    

    while(1)
      {
	const unsigned int slots = atomicBuildRecordStack.usedSlots();

	if (slots >= NUM_TOTAL_CORES) break;

	if (slots >= MAX_BUILD_RECORD_STACK_ENTRIES) break;
	ParallelQBVHBuilder::pushLargestElementToTop();

	if (slots)
	  {
	    const unsigned int sizeTopElement = atomicBuildRecordStack.t[slots-1].items();
	    if ( unlikely(sizeTopElement < BUILD_RECORD_SPLIT_THRESHOLD))
	      break;
	    if ( unlikely(slots >= ParallelQBVHBuilder::numActiveLocalWorkQueues &&
			  sizeTopElement < BUILD_RECORD_SPLIT_THRESHOLD*2) )
	      break;
	  }

	if (atomicBuildRecordStack.pop_nolock(record) == false) break;
	unsigned int localNodeID = 0;
	unsigned int localNodeIDs = 0;
	const unsigned int numSplits = recursePtr(record,PUSH_GLOBAL,0,localNodeID,localNodeIDs);
      }

    TIMER(cycles = (unsigned long)timer.stop());
    TIMER(std::cout << "GENERATE " << atomicBuildRecordStack.usedSlots() << " QBVH NODES " << cycles / 1000000.0f << "MCYCLES " << std::endl);
    }


    // ===============================================================
    CHECK(
	  std::cout << "CHECK BUILD RECORDS 1" << std::endl;
	  for (int i=0;i<atomicBuildRecordStack.usedSlots();i++) checkBuildRecord(atomicBuildRecordStack.get(i));
	  std::cout << "DONE" << std::endl;
	  );


    // ===============================================================

    TIMER(timer.start());  

    quicksort_decending<BuildRecord>(atomicBuildRecordStack.get_addr(0),0,atomicBuildRecordStack.usedSlots()-1);
    transferBuildRecordsFromGlobalToLocalStack();

    TIMER(cycles = (unsigned long)timer.stop());
    TIMER(std::cout << "sort and transfer build records " << cycles / 1000000.0f << " Mcycles " << std::endl);
  
    CHECK(
	  DBG_PRINT(atomicBuildRecordStack.usedSlots());
	  for (int i=0;i<atomicBuildRecordStack.usedSlots();i++) DBG_PRINT(atomicBuildRecordStack.t[i].items())
	  );
 
    // ===============================================================

    TIMER(timer.start());  
    QBVHTaskScheduler::dispatchTask( ParallelQBVHBuilder::thread_build_local );
    TIMER(cycles = (unsigned long)timer.stop());
    TIMER(std::cout << "thread_build " << cycles << " cycles " << std::endl);

    // ===============================================================

    if (atomicID.val() % 4 != 0) FATAL("number of nodes % 4 != 0");

    qbvh_nodes = atomicID.val() >> 2;

    // ===============================================================

    TIMER(timer.start());
    QBVHTaskScheduler::dispatchTask(ParallelQBVHBuilder::thread_createTriangleAccel);
    TIMER(cycles = (unsigned long)timer.stop());
    TIMER(std::cout << "Created triangle acceleration in " << cycles << " cycles " << cycles/1000000.0f << " mcycles" << std::endl << std::flush);

    CHECK( DBG_PRINT(qbvh_nodes) );
    CHECK( std::cout << "SAH BUILDER DONE" << std::endl << std::flush );
  }

  // =================================================================
  // =================================================================
  // =================================================================

void ParallelQBVHBuilder::checkBuildRecord(const BuildRecord &current)
  {
    AABB check_box;
    AABB box;
    check_box.setEmpty();
    box.setEmpty();
    current.cs_AABB.storeSceneAABB((float*)&box);

    for (unsigned int i=current.begin;i<current.end;i++) 
      {
	check_box.extend(aabb[i]);
	if (!box.enclose(aabb[i]))
	  {
	    DBG_PRINT(current);
	    DBG_PRINT(i);
	    DBG_PRINT(aabb[i]);
	    FATAL("check build record");
	  }
      }
    if (!equal(check_box,box)) 
      {
	DBG_PRINT(current);
	DBG_PRINT(check_box);
	DBG_PRINT(box);
	FATAL("check build record");
      }
  }

  bool ParallelQBVHBuilder::checkRange(const unsigned int begin,
				       const unsigned int end,
				       const unsigned int side)
  {
    const unsigned int bestSplit    = splitData.bestSplit;
    const unsigned int bestSplitDim = splitData.bestSplitDim;
    const mic_f centroidBoundsMin_2 = splitData.centroidBoundsMin_2;
    const mic_f scale               = splitData.scale;

    MIC_ALIGN Centroid_Scene_AABB box;
    box.reset();

    if (side == 0)
      {
	for (int i=begin;i<end;i++)
	  if (getBinSlot(aabb+i,bestSplitDim,centroidBoundsMin_2,scale) >= bestSplit)
	    {
	      DBG_PRINT(begin);
	      DBG_PRINT(end);
	      DBG_PRINT(side);
	      std::cout << " ERROR " << i << " " << getBinSlot(aabb+i,bestSplitDim,centroidBoundsMin_2,scale) << " " << bestSplit << std::endl;
	      return false;
	    }
	  else
	    box.extend(aabb[i]);
      }
    else
      {
	for (int i=begin;i<end;i++)
	  if (getBinSlot(aabb+i,bestSplitDim,centroidBoundsMin_2,scale) < bestSplit)
	    {
	      DBG_PRINT(begin);
	      DBG_PRINT(end);
	      DBG_PRINT(side);
	      std::cout << " ERROR " << i << " " << getBinSlot(aabb+i,bestSplitDim,centroidBoundsMin_2,scale) << " " << bestSplit << std::endl;
	      return false;
	    }
	  else
	    box.extend(aabb[i]);
      }
    return true;
  }

};
