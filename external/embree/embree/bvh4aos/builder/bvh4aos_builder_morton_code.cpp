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

#include "bvh4aos_globals.h"
#include "bvh4aos_builder_util.h"
#include "bvh4aos_builder_common.h"

#define CHECK(x) 
#define DBG(x) 
#define DBG_CHECK(x) 
#define TIMER(x) 

#define L1_PREFETCH_ITEMS 8
#define L2_PREFETCH_ITEMS 44

#define LATTICE_BITS_PER_DIM 10
#define LATTICE_SIZE_PER_DIM ((unsigned int)1 << LATTICE_BITS_PER_DIM)

#define MAX_TOP_LEVEL_BINS 1024

#define TOP_LEVEL_SPLIT_POS 10
#define MORTON_IDS_PER_BLOCK 8
#define MAX_REFIT_SUBTREES 512

#define RADIX_BITS    8
#define RADIX_BUCKETS (1 << RADIX_BITS)
#define RADIX_BUCKETS_MASK (RADIX_BUCKETS-1)

#define LOCAL_BUILD_RECORD_PUSH_THRESHOLD 128
#define REFIT_TRIANGLES_PER_THREAD 64

namespace embree
{


  MortonID32Bit *ParallelQBVHBuilderMortonCode::mortonID[2]       = { NULL, NULL };
  RadixNode     *ParallelQBVHBuilderMortonCode::radixTree         = NULL;
  unsigned int   ParallelQBVHBuilderMortonCode::numMortonIDBlocks = 0;

  static unsigned int subTrees  = 0;
  static unsigned int changes01 = 0;

  MIC_ALIGN static unsigned int subTreeIDs[MAX_REFIT_SUBTREES];
  MIC_ALIGN static unsigned int   nodeIDs[MAX_TOP_LEVEL_BINS];
  MIC_ALIGN static MortonID32Bit change01[MAX_TOP_LEVEL_BINS];
  MIC_ALIGN static unsigned int matrix[MAX_MIC_THREADS][RADIX_BUCKETS];
  MIC_ALIGN static unsigned int sum_matrix[MAX_MIC_THREADS];



  _INLINE std::string format_binary(unsigned int x)
  {
    static char b[33];
    b[32] = '\0';

    for (int z = 0; z < 32; z++) {
      b[31-z] = ((x>>z) & 0x1) ? '1' : '0';
    }

    return b;
  }


  // ==================================================================================================================

  void ParallelQBVHBuilderMortonCode::thread_radixsort_blocked(const unsigned int threadID)
  {
    const unsigned int WORKERS = NUM_TOTAL_THREADS;
    const unsigned int workerID = threadID;
    const unsigned int primitives_per_block = splitData.work_items_per_thread;

    const unsigned int startBlockID = (workerID * primitives_per_block);
    const unsigned int endBlockID   = (workerID != WORKERS-1) ? (startBlockID + primitives_per_block) : numMortonIDBlocks; 
    const unsigned int startID      = startBlockID * MORTON_IDS_PER_BLOCK;
    const unsigned int endID        = endBlockID   * MORTON_IDS_PER_BLOCK;

    MIC_ALIGN unsigned int offset[RADIX_BUCKETS];

    for (unsigned long b=0;b<4;b++)
      {
#pragma unroll(RADIX_BUCKETS/16)
	for (unsigned long j=0;j<RADIX_BUCKETS;j+=16)
	  prefetch<PFHINT_L2EX>(&matrix[threadID][j]);

	syncThreads(workerID,WORKERS);

	const MortonID32Bit *__restrict source = (MortonID32Bit *)&mortonID[(b%2)][0];
	MortonID32Bit *__restrict dest   = (MortonID32Bit *)&mortonID[1-(b%2)][0];

#pragma unroll(RADIX_BUCKETS/16)
	for (unsigned long j=0;j<RADIX_BUCKETS;j+=16)
	  prefetch<PFHINT_L1EX>(&matrix[threadID][j]);

	for (unsigned long i=0;i<RADIX_BUCKETS;i++)
	  matrix[threadID][i] = 0;

	// ----------------------------------------------------------

	for (unsigned long i=startID;i<endID;i+=MORTON_IDS_PER_BLOCK)
	  {
	    prefetch<PFHINT_L1>(&source[i+L1_PREFETCH_ITEMS]);
	    prefetch<PFHINT_L2>(&source[i+L2_PREFETCH_ITEMS]);


#pragma unroll(MORTON_IDS_PER_BLOCK)
	    for (unsigned long j=0;j<MORTON_IDS_PER_BLOCK;j++)
	      {
		const unsigned long index = source[i+j].getByte(b);
		assert(index < RADIX_BUCKETS);
		matrix[threadID][index]++;
	      }

	  }

	syncThreads(workerID,WORKERS);

      

	if (likely(NUM_TOTAL_THREADS >=64))
	  {
	    const unsigned long coreID = workerID >> 2;
	    if (coreID < 16 && (threadID % 4 == 0))
	      {
		mic_i column = mic_i::zero();      

		for (unsigned long i=0;i<NUM_TOTAL_THREADS;i++)
		  prefetch<PFHINT_L2EX>(&matrix[i][coreID * 16]);

		for (unsigned long i=0;i<NUM_TOTAL_THREADS-1;i++)
		  {
		    prefetch<PFHINT_L1EX>(&matrix[i+1][coreID * 16]);
		    const mic_i old_column = column;
		    column += load16i((const int *)&matrix[i][coreID * 16]);
		    store16i(&matrix[i][coreID * 16],old_column);	      
		    evictL1(&matrix[i][coreID * 16]);
		  }
		{
		  const mic_i old_column = column;
		  column += load16i((const int *)&matrix[NUM_TOTAL_THREADS-1][coreID * 16]);
		  store16i(&matrix[NUM_TOTAL_THREADS-1][coreID * 16],old_column);
		  evictL1(&matrix[NUM_TOTAL_THREADS-1][coreID * 16]);
		}

		store16i_ngo(&sum_matrix[coreID * 16],column);
	      }
	  }
	else
	  {
	    if (threadID == 0)
	      for (unsigned int coreID=0;coreID<16;coreID++)
		{
		  mic_i column = mic_i::zero();      

		  for (unsigned long i=0;i<NUM_TOTAL_THREADS;i++)
		    prefetch<PFHINT_L2EX>(&matrix[i][coreID * 16]);

		  for (unsigned long i=0;i<NUM_TOTAL_THREADS-1;i++)
		    {
		      prefetch<PFHINT_L1EX>(&matrix[i+1][coreID * 16]);
		      const mic_i old_column = column;
		      column += load16i((const int *)&matrix[i][coreID * 16]);
		      store16i(&matrix[i][coreID * 16],old_column);	      
		      evictL1(&matrix[i][coreID * 16]);
		    }
		  {
		    const mic_i old_column = column;
		    column += load16i((const int *)&matrix[NUM_TOTAL_THREADS-1][coreID * 16]);
		    store16i(&matrix[NUM_TOTAL_THREADS-1][coreID * 16],old_column);
		    evictL1(&matrix[NUM_TOTAL_THREADS-1][coreID * 16]);
		  }

		  store16i_ngo(&sum_matrix[coreID * 16],column);
		}
	  }
	    

	for (unsigned long i=0;i<RADIX_BUCKETS;i+=16)
	  prefetch<PFHINT_L1EX>(&offset[i]);

	syncThreads(workerID,WORKERS);

#pragma unroll(16)
	for (unsigned long i=0;i<RADIX_BUCKETS;i+=16)
	  prefetch<PFHINT_L1>(&sum_matrix[i]);
      
	offset[0] = 0;
	for (unsigned long i=1;i<RADIX_BUCKETS;i++)    
	  offset[i] = offset[i-1] + sum_matrix[i-1];

      
#pragma unroll(RADIX_BUCKETS)
	for (unsigned long j=0;j<RADIX_BUCKETS;j++)
	  {
	    offset[j] += matrix[threadID][j];
	  }
      
	for (unsigned long i=startID;i<endID;i+=MORTON_IDS_PER_BLOCK)
	  {
	    prefetch<PFHINT_L1>(&source[i+L1_PREFETCH_ITEMS]);
	    prefetch<PFHINT_L2>(&source[i+L2_PREFETCH_ITEMS]);

#pragma unroll(MORTON_IDS_PER_BLOCK)
	    for (unsigned long j=0;j<MORTON_IDS_PER_BLOCK;j++)
	      {	  
		const unsigned long index = source[i+j].getByte(b);
		assert(index < RADIX_BUCKETS);
		dest[offset[index]] = source[i+j];
		prefetch<PFHINT_L2EX>(&dest[offset[index]+L1_PREFETCH_ITEMS]);
		offset[index]++;
	      }
	    evictL2(&source[i]);
	  }
      }
  }

  // ==================================================================================================================
  // ==================================================================================================================
  // ==================================================================================================================

  _INLINE unsigned long delta(const MortonID32Bit *__restrict__ const m,
			      const unsigned long begin,
			      const unsigned long end)
  {
    const unsigned int code_start = m[begin].code;
    const unsigned int code_end   = m[end].code;
    const unsigned int bitpos     = clz(code_start^code_end);
    return bitpos;
  }



  _INLINE long getSplitPosition(const MortonID32Bit *__restrict__ const m,
				const unsigned long begin,
				const unsigned long end)
  {
    long index = -1;

    const unsigned int code_start = m[begin].code;
    const unsigned int code_end   = m[end-1].code;
    const unsigned int bitpos     = clz(code_start^code_end);
    //assert(bitpos == bitpos2);

    if (bitpos > TOP_LEVEL_SPLIT_POS) return index;

    if (unlikely(bitpos == 32)) {
      index = begin + (end-begin+1)/2; // all have the same morton code
      std::cout << "WARNING" << std::endl;
    }
    else
      {
	const unsigned int bitpos_diff = 31-bitpos;
	const unsigned int bitmask = (unsigned int)1 << (bitpos_diff);

	unsigned long scan_start = begin;
	unsigned long scan_end   = end;
	while(1)
	  {
	    assert(scan_start < scan_end);
	    if (unlikely(scan_start + 1 == scan_end))
	      {
		assert( (m[scan_start+0].code & bitmask) == 0 );
		assert( (m[scan_start+1].code & bitmask) == bitmask );
		break;
	      }

	    const unsigned long mid = ((scan_end+scan_start) >> 1);
	    if ((m[mid].code & bitmask) == 0)
	      {
		scan_start = mid;
	      
	      }
	    else
	      {
		scan_end = mid;
	      }
	  }
	index = scan_start + 1;

#if defined(DEBUG)
	assert(index != -1);
	for (unsigned long i=index;i<end;i++)
	  assert((m[i].code & bitmask) == bitmask);
#endif
      }
    assert(index != -1);
    return index;
  }

				
  // ==================================================================================================================
  // ==================================================================================================================
  // ==================================================================================================================

  _INLINE void split2_morton(BuildRecord &current,
			     BuildRecord &left,
			     BuildRecord &right,
			     const MortonID32Bit *__restrict__ const m)
  {
    if (unlikely(current.items() <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD))
      {
	current.createLeaf();
	return;
      }

    const unsigned int code_start = m[current.begin].code;
    const unsigned int code_end   = m[current.end-1].code;
    const unsigned int bitpos = clz(code_start^code_end);

    long index = -1;

    if (unlikely(bitpos == 32)) 
      {
	index = current.begin + current.items()/2; // all have the same morton code
      }
    else
      {
	const unsigned int bitpos_diff = 31-bitpos;
	const unsigned int bitmask = (unsigned int)1 << (bitpos_diff);

	unsigned long start = current.begin;
	unsigned long end   = current.end;
	while(1)
	  {
	    const unsigned long mid = ((end+start) >> 1);
	    prefetch<PFHINT_L1>(&m[mid]);
	    assert(start < end);

	    if (unlikely(start + 1 == end))
	      {
		assert( (m[start+0].code & bitmask) == 0 );
		assert( (m[start+1].code & bitmask) == bitmask );
		break;
	      }

	    if ((m[mid].code & bitmask) == 0)
	      {
		start = mid;
	      
	      }
	    else
	      {
		end = mid;
	      }
	  }
	index = start + 1;

#if defined(DEBUG)      
	assert(index != -1);
	for (unsigned long i=index;i<current.end;i++)
	  assert((m[i].code & bitmask) == bitmask);
#endif

      }
    left.init(current.begin,index,0);
    right.init(index,current.end,0);
  }


  // ==================================================================================================================
  // ==================================================================================================================
  // ==================================================================================================================

  unsigned int ParallelQBVHBuilderMortonCode::recurseMortonID(BuildRecord &current, 
							      const unsigned int mode, 
							      const unsigned int threadID,
							      unsigned int &localNodeID,
							      unsigned int &localNodeIDs)
  {
    MIC_ALIGN BuildRecord record[4];
    MIC_ALIGN BuildRecord left,right;

    unsigned int numSplits = 1;

    const MortonID32Bit *__restrict__ const m = mortonID[0];

    prefetch<PFHINT_L1>(&m[current.begin]);

    if (current.items() <= QBVH_BUILDER_LEAF_ITEM_THRESHOLD) 
      {
      }
    else
      {
	record[0] = current;

	left.prefetchL1Ex();
	right.prefetchL1Ex();

	while(1)
	  {
	    bool couldSplit = false;
	    if (numSplits < 4)
	      {
		int index = -1;
		unsigned int items = 0;
		for (unsigned int i=0;i<numSplits;i++)
		  {
		    if (!record[i].isLeaf())
		      {
			if (record[i].items() >= items)
			  {
			    items = record[i].items();
			    index = i;
			  }
		      }
		  }
		assert(index < 4);
		if (index != -1)
		  {
		    split2_morton(record[index],left,right,m);

		    assert(numSplits > 0);
		    record[numSplits].prefetchL1Ex();

		    if ( !record[index].isLeaf() )
		      {
			record[index] = record[numSplits-1];
			numSplits--;

			record[numSplits+0] = left;
			record[numSplits+1] = right;
			numSplits+=2;

			prefetch<PFHINT_L2>(globalBuildTrianglePtr + m[left.begin].index);
			prefetch<PFHINT_L2>(globalBuildTrianglePtr + m[right.begin].index);

			couldSplit = true;
		      }
		    else 
		      {
			continue;
		      }
		  }
	      }
	    // ==================================================
	    if (!couldSplit) break;
	  }
      }

    // could not generate any split, propagate leaf back to parent
    if ( numSplits == 1 )
      {
	if(current.items() <= 4)
	  {
	    DBG(std::cout << "rollback forced leaf with " << current.items() << " items" << std::endl);
	    DBG(std::cout << "couldn't split -> propagate back to parent" << std::endl);

#if defined(DEBUG)
	    if (unlikely(!(mode == RECURSE || mode == PUSH_LOCAL))) std::cout << "warning leaf" << std::endl;
#endif

	    current.createLeaf();

	    const unsigned int items  = current.items();
	    const unsigned int offset = current.begin;
	    assert(items  != 0);


	    getLeafBoundsAndComputeAccelerationData(offset,items,current.nodeID);
	    globalNodePtr[current.nodeID].createLeaf(offset,items,0);

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
      FATAL("not enough nodes allocated");

    prefetch<PFHINT_L1EX>(&globalNodePtr[currentIndex+0]);
    prefetch<PFHINT_L1EX>(&globalNodePtr[currentIndex+2]);
    prefetch<PFHINT_L2EX>(&globalNodePtr[currentIndex+4]);
    prefetch<PFHINT_L2EX>(&globalNodePtr[currentIndex+6]);

    current.childrenID = currentIndex;

    // ==== init all four nodes ====
    const mic_f init_node = upconv16f((float*)ParallelQBVHBuilder::initQBVHNode);


    for (int i=0;i<4;i++)
      if (record[i].items() < QBVH_BUILDER_LEAF_ITEM_THRESHOLD)
	record[i].createLeaf();

    unsigned int nodes_to_process = 0;
    unsigned int node_to_process[4];

    const unsigned int coreID = threadID >> 2;

    for (unsigned int i=0;i<numSplits;i++)
      {
	assert(record[i].items());

	if (!record[i].isLeaf())
	  {
	    record[i].nodeID = currentIndex+i;

	    if (record[i].items() > LOCAL_BUILD_RECORD_PUSH_THRESHOLD &&
		localAtomicBuildRecordStack[coreID].usedSlots() < NUM_LOCAL_ATOMIC_BUILD_RECORDS-1)
	      {		  
		while (localAtomicBuildRecordStack[coreID].push(record[i]) == false)
		  {
		    delayThread(1024);
		    DBG(std::cout << "LOCAL ENQUEUE ERROR " << coreID << std::endl);
		  }
	      }
	    else
	      node_to_process[nodes_to_process++] = i;
	  }
	else
	  {
	    const unsigned int items  = record[i].items();
	    const unsigned int offset = record[i].begin;
	  
	    assert (items <= 4);
	    getLeafBoundsAndComputeAccelerationData(offset,items,currentIndex+i);

	    DBG(std::cout << "create leaf to node " << currentIndex+i << " items " << items << " offset " << offset << std::endl);

	    globalNodePtr[currentIndex+i].createLeaf(offset,items,0);
	  }
      
      }

    for (unsigned int i=numSplits;i<4;i++)
      globalNodePtr[currentIndex+i] = *(BVHNode*)&initQBVHNode[0];
      

    globalNodePtr[current.nodeID].createNode(current.childrenID,numSplits,current.items());

    assert(current.childrenID > 0);
    assert(numSplits > 0);

    if (nodes_to_process > 0)
      for (unsigned int i=0;i<nodes_to_process;i++)
	{
	  const unsigned int index = node_to_process[i];
	  recurseMortonID(record[index],mode,threadID,localNodeID,localNodeIDs);
	}

    return numSplits;  
  }

  // ==================================================================================================================
  // ==================================================================================================================
  // ==================================================================================================================

  void ParallelQBVHBuilderMortonCode::thread_getSceneBounds(const unsigned int threadID)
  {
    const unsigned int WORKERS = NUM_TOTAL_THREADS;
    const unsigned int workerID = threadID;
    unsigned int primitives_per_block = splitData.work_items_per_thread;

    const unsigned long startID = workerID * primitives_per_block;
    const unsigned long endID  = (workerID != WORKERS-1) ? (startID + primitives_per_block) : vertices;

    Centroid_Scene_AABB cs_AABB;
    cs_AABB.reset();
    mic_f centroidMinAABB = cs_AABB.minCentroidAABB();
    mic_f centroidMaxAABB = cs_AABB.maxCentroidAABB();


    const Vec3fa * __restrict__ const vptr = globalVertexPtr;

    for (unsigned long i=startID;i<endID;i++)
      {
	prefetch<PFHINT_NT>(&vptr[i+4]);
	prefetch<PFHINT_L2>(&vptr[i+L2_PREFETCH_ITEMS]);
	const mic_f v = upconv4f((float*)&vptr[i]); 
	centroidMinAABB = _min(centroidMinAABB,v);
	centroidMaxAABB = _max(centroidMaxAABB,v);
      }

    cs_AABB.setMinCentroidAABB(centroidMinAABB);
    cs_AABB.setMaxCentroidAABB(centroidMaxAABB);

    global_cs_AABB.extend_atomic_centroid_bounds(cs_AABB);
  }


  // ==================================================================================================================
  // ==================================================================================================================
  // ==================================================================================================================



  MIC_ALIGN static int mIndexTable[16] = { 0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7 };

  void ParallelQBVHBuilderMortonCode::thread_computeMortonIDs(const unsigned int threadID)
  {
    MIC_ALIGN mic_f centroid_x;
    MIC_ALIGN mic_f centroid_y;
    MIC_ALIGN mic_f centroid_z;

    prefetch<PFHINT_L1EX>(&centroid_x);
    prefetch<PFHINT_L1EX>(&centroid_y);
    prefetch<PFHINT_L1EX>(&centroid_z);

    const unsigned int WORKERS = NUM_TOTAL_THREADS;
    const unsigned int workerID = threadID;
    const unsigned int primitives_per_block = splitData.work_items_per_thread;

    
    const unsigned int startBlockID = (workerID * primitives_per_block);
    const unsigned int endBlockID   = (workerID != WORKERS-1) ? (startBlockID + primitives_per_block) : numMortonIDBlocks-1; 
    const unsigned int startID      = startBlockID * MORTON_IDS_PER_BLOCK;
    const unsigned int endID        = endBlockID * MORTON_IDS_PER_BLOCK; 
  
    assert(startID % 8 == 0);
    assert(endID % 8   == 0);

    const Vec3fa          * __restrict__ const vptr = globalVertexPtr;
    const BuilderTriangle * __restrict__       tptr = globalBuildTrianglePtr;

    const mic_f boundsMin = global_cs_AABB.minCentroidAABB();
    const mic_f diagonal  = global_cs_AABB.centroidDiagonal();
    const mic_f scale     = sel(nz(diagonal),rcp(diagonal) * mic_f(LATTICE_SIZE_PER_DIM * 0.99f),mic_f::zero());

    MortonID32Bit *__restrict__ const m = mortonID[0];


    const mic_f bMin_x(swAAAA(boundsMin));
    const mic_f bMin_y(swBBBB(boundsMin));
    const mic_f bMin_z(swCCCC(boundsMin));

    const mic_f dScale_x(swAAAA(scale));
    const mic_f dScale_y(swBBBB(scale));
    const mic_f dScale_z(swCCCC(scale));

    centroid_x = mic_f::zero();
    centroid_y = mic_f::zero();
    centroid_z = mic_f::zero();


    const mic_i iTable = load16i(mIndexTable);

    for (unsigned long i=startID;i<endID;i+=MORTON_IDS_PER_BLOCK)
      {
	#pragma unroll(MORTON_IDS_PER_BLOCK)
	for (unsigned long j=0;j<MORTON_IDS_PER_BLOCK;j++)
	  {
	    const mic_f v0 = upconv4f((float*)&vptr[tptr[i+j].v0]); 
	    const mic_f v1 = upconv4f((float*)&vptr[tptr[i+j].v1]);
	    const mic_f v2 = upconv4f((float*)&vptr[tptr[i+j].v2]);
	  
	    const mic_f b_min = _min(v0,_min(v1,v2));
	    const mic_f b_max = _max(v0,_max(v1,v2));
	    const mic_f c     = (b_min+b_max) * mic_f(0.5f); 
	    compactustore16f_low(0x1,&centroid_x[2*j+0],c);
	    compactustore16f_low(0x2,&centroid_y[2*j+0],c);
	    compactustore16f_low(0x4,&centroid_z[2*j+0],c);
	  }
	evictL2(&tptr[i+0]);
	const mic_i index = iTable + i;
	evictL2(&tptr[i+4]);

	const mic_i binID_x = mic_i((centroid_x - bMin_x)*dScale_x);
	const mic_i binID_y = mic_i((centroid_y - bMin_y)*dScale_y);
	const mic_i binID_z = mic_i((centroid_z - bMin_z)*dScale_z);

	prefetch<PFHINT_NT>(&tptr[i+8]);
	prefetch<PFHINT_NT>(&tptr[i+16]);

	const mic_i code  = bitInterleave(binID_x,binID_y,binID_z);
	const mic_i final = sel(0x5555,code,index);
	store16i_ngo(&m[i],final);
      }

    if (workerID == WORKERS-1)
      {
	const unsigned int last_startID = (numMortonIDBlocks-1) * MORTON_IDS_PER_BLOCK;

	centroid_x = mic_f::zero();
	centroid_y = mic_f::zero();
	centroid_z = mic_f::zero();
	unsigned long j = 0;
	for (unsigned long i=last_startID;i<triangles;i++)
	  {
	    const mic_f v0 = upconv4f((float*)&vptr[tptr[i].v0]); 
	    const mic_f v1 = upconv4f((float*)&vptr[tptr[i].v1]);
	    const mic_f v2 = upconv4f((float*)&vptr[tptr[i].v2]);
	  
	    const mic_f b_min = _min(v0,_min(v1,v2));
	    const mic_f b_max = _max(v0,_max(v1,v2));
	    const mic_f c     = (b_min+b_max) * mic_f(0.5f); 
	    compactustore16f_low(0x1,&centroid_x[2*j+0],c);
	    compactustore16f_low(0x2,&centroid_y[2*j+0],c);
	    compactustore16f_low(0x4,&centroid_z[2*j+0],c);
	    j++;
	  }
	const mic_i index = iTable + last_startID;
	const mic_i binID_x = mic_i((centroid_x - bMin_x)*dScale_x);
	const mic_i binID_y = mic_i((centroid_y - bMin_y)*dScale_y);
	const mic_i binID_z = mic_i((centroid_z - bMin_z)*dScale_z);
	const mic_i code  = bitInterleave(binID_x,binID_y,binID_z);
	const mic_i final = sel(0x5555,code,index);
	store16i_ngo(&m[last_startID],final);
      }

  }

  // ==================================================================================================================
  // ==================================================================================================================
  // ==================================================================================================================


  void ParallelQBVHBuilderMortonCode::thread_getTaskIntervals(const unsigned int threadID)
  {
    const unsigned int WORKERS = NUM_TOTAL_THREADS;
    const unsigned int workerID = threadID;
    const unsigned int primitives_per_block = splitData.work_items_per_thread;

    const unsigned int startBlockID = (workerID * primitives_per_block);
    const unsigned int endBlockID   = (workerID != WORKERS-1) ? (startBlockID + primitives_per_block) : numMortonIDBlocks; 
    const unsigned int startID      = startBlockID * MORTON_IDS_PER_BLOCK;
    const unsigned int endID        = endBlockID   * MORTON_IDS_PER_BLOCK;

    const MortonID32Bit *__restrict__ const m = mortonID[0];

    for (unsigned long i=startID;i<endID;i+=MORTON_IDS_PER_BLOCK)
      {
	const mic_i code_start =  load16i((int*)&m[i+0]);
	const mic_i code_end   = uload16i((int*)&m[i+1]);
	const mic_i diff       = bitReverse(code_start^code_end);

	prefetch<PFHINT_NT>(&m[i+1*MORTON_IDS_PER_BLOCK]);
	prefetch<PFHINT_NT>(&m[i+2*MORTON_IDS_PER_BLOCK]);

	for (unsigned long j=0;j<MORTON_IDS_PER_BLOCK;j++)
	  if (unlikely(bsf32(diff[2*j]) <= TOP_LEVEL_SPLIT_POS))
	    {
	      if (unlikely(i+j+1 >= triangles)) continue;
	      assert(delta(m,i+j,i+j+1) <= TOP_LEVEL_SPLIT_POS);
	      assert(i+j < triangles);
	      const unsigned int index = ParallelQBVHBuilder::atomicID.inc();
	      change01[index].code  = m[i+j+1].code;
	      change01[index].index = i+j+1;
	    }	  
      }

    if (threadID == 0)
      {
	for (unsigned int i=0;i<NUM_TOTAL_CORES;i++)
	  {
	    for (unsigned int j=0;j<4;j++)
	      localAtomicBuildRecordStack[i].get(j).prefetchL2Ex();
	    localAtomicBuildRecordStack[i].init();
	  }
      }
    
  }

  // ==================================================================================================================
  // ==================================================================================================================
  // ==================================================================================================================

  void ParallelQBVHBuilderMortonCode::thread_createBuildRecordsFromTaskIntervals(const unsigned int threadID)
  {
    MIC_ALIGN BuildRecord rec;

    const unsigned int numIntervals = changes01-1;
    const unsigned int coreID       = threadID >> 2;
  
    unsigned int i = (threadID >> 2) + (threadID % 4) * NUM_TOTAL_CORES;
    while(i < numIntervals)
      {
	const unsigned int start = change01[i].index;
	const unsigned int end   = change01[i+1].index;
	rec.init(global_cs_AABB,start,end);
	rec.nodeID = nodeIDs[i];
	localAtomicBuildRecordStack[coreID].push(rec);
	i += NUM_TOTAL_THREADS;
      }
  }

  void ParallelQBVHBuilderMortonCode::createBuildRecordsFromTaskIntervals()
  {
    MIC_ALIGN unsigned int distribution[64];
    MIC_ALIGN unsigned int entries[64];

    for (unsigned int i=0;i<numActiveLocalWorkQueues;i++) 
      {
	distribution[i] = 0;
	entries[i] = 0;
      }

    for (unsigned int i=0;i<changes01-1;i++)
      {
	MIC_ALIGN BuildRecord rec;
	const unsigned int start = change01[i].index;
	const unsigned int end   = change01[i+1].index;


	rec.init(global_cs_AABB,start,end,0);
	rec.nodeID = nodeIDs[i];

#if defined(DEBUG)
	for (unsigned long j=start+1;j<end;j++)
	  {
	    if (mortonID[0][j] < mortonID[0][j-1] || delta(mortonID[0],j,j-1) < TOP_LEVEL_SPLIT_POS)
	      {
		std::cout << "ERROR" << std::endl;
		DBG_PRINT(j);
		DBG_PRINT(j-1);

		DBG_PRINT(delta(mortonID[0],j,j-1) );
		DBG_PRINT( mortonID[0][j] );
		DBG_PRINT( mortonID[0][j-1] );
		FATAL("HERE");
	      }
	  }
#endif


	unsigned int smallest_index = 0;

	for (unsigned int j=0;j<numActiveLocalWorkQueues;j++)
	  if (entries[j] < NUM_LOCAL_ATOMIC_BUILD_RECORDS)
	    {
	      smallest_index = j;
	      break;
	    }

	unsigned int smallest = distribution[smallest_index];

	for (unsigned int j=smallest_index+1;j<numActiveLocalWorkQueues;j++)
	  if (distribution[j] < smallest && entries[j] < NUM_LOCAL_ATOMIC_BUILD_RECORDS)
	    {
	      smallest = distribution[j];
	      smallest_index = j;
	    }

	distribution[smallest_index] += rec.items();
	entries[smallest_index]++;

	localAtomicBuildRecordStack[smallest_index].push_nolock(rec);
      }
#if 0
    for (int i=0;i<numActiveLocalWorkQueues;i++)
      {
	unsigned int total = 0;
	for (int j=0;j<localAtomicBuildRecordStack[i].usedSlots();j++)
	  total += localAtomicBuildRecordStack[i].t[j].items();
	std::cout << "core " << i << " -> " << total << std::endl;
      }
#endif

  }

  // ==================================================================================================================
  // ==================================================================================================================
  // ==================================================================================================================

  void ParallelQBVHBuilderMortonCode::createTopLevel(const unsigned long start,
						     const unsigned long end,
						     const unsigned long fatherID,
						     const MortonID32Bit *__restrict__ const m,
						     unsigned long &i_index)
  {
    const unsigned int items = end - start;
    assert(items>0);
    const mic_f init_node = upconv16f((float*)ParallelQBVHBuilder::initQBVHNode);

  
    if (items == 1) 
      {
	nodeIDs[i_index++] = fatherID;
      }
    else
      {

	const unsigned int ID = ParallelQBVHBuilder::atomicID.add(4);

	store16f_ngo(&globalNodePtr[ID+0],init_node);
	store16f_ngo(&globalNodePtr[ID+2],init_node);

	unsigned int children = 0;

	const long splitm = getSplitPosition(m,start,end);

	assert(splitm != -1);
	assert(splitm < end);


	const long splitl = getSplitPosition(m,start,splitm);
	if (splitl != -1)
	  {	  
	    assert(splitl > start);
	    assert(splitl < splitm);
	    createTopLevel(start, splitl,ID+children,m,i_index);
	    children++;
	    createTopLevel(splitl,splitm,ID+children,m,i_index);
	    children++;
	  }
	else
	  {
	    createTopLevel(start, splitm,ID+children,m,i_index);
	    children++;
	  }

	const long splitr = getSplitPosition(m,splitm,end);
	if (splitr != -1)
	  {      
	    assert(splitr > splitm);
	    assert(splitr < end);
	    createTopLevel(splitm,splitr,ID+children,m,i_index);
	    children++;
	    createTopLevel(splitr,   end,ID+children,m,i_index);
	    children++;
	  }
	else
	  {
	    createTopLevel(splitm,end,ID+children,m,i_index);	        
	    children++;
	  }

	assert(children >=2);
	assert(children <=4);
	globalNodePtr[fatherID].createNode(ID,children,triangles);
      }
  }


  // ==================================================================================================================
  // ==================================================================================================================
  // ==================================================================================================================

  void ParallelQBVHBuilderMortonCode::handleLargeLeaves(const unsigned int nodeID,
							const unsigned int start,
							const unsigned int items)
  {
    assert(items > 4);
    const unsigned int ID = ParallelQBVHBuilder::atomicID.add(4);
      
    for (unsigned int i=0;i<4;i++)
      globalNodePtr[ID+i] = *(BVHNode*)&initQBVHNode[0];

    const unsigned int blocks = (items+3)/4;

    globalNodePtr[nodeID].createNode(ID,min((int)blocks,4),items);


    getLeafBoundsAndComputeAccelerationData(start,4,ID+0);
    globalNodePtr[ID+0].createLeaf(start,4.0);

    if (items <= 8)
      {
	getLeafBoundsAndComputeAccelerationData(start+4,items-4,ID+1);
	globalNodePtr[ID+1].createLeaf(start+4,items-4.0);	  
      }
    else if (items <= 12)
      {
	getLeafBoundsAndComputeAccelerationData(start+4,4,ID+1);
	globalNodePtr[ID+1].createLeaf(start+4,4.0);	  
	getLeafBoundsAndComputeAccelerationData(start+8,items-8,ID+2);
	globalNodePtr[ID+2].createLeaf(start+8,items-8.0);	  
      }
    else if (items <= 16)
      {
	getLeafBoundsAndComputeAccelerationData(start+4,4,ID+1);
	globalNodePtr[ID+1].createLeaf(start+ 4,4.0);	  
	getLeafBoundsAndComputeAccelerationData(start+8,4,ID+2);
	globalNodePtr[ID+2].createLeaf(start+ 8,4.0);	  
	getLeafBoundsAndComputeAccelerationData(start+12,items-12,ID+3);
	globalNodePtr[ID+3].createLeaf(start+12,items-12.0);	  
      }
    else
      {
	getLeafBoundsAndComputeAccelerationData(start+4,4,ID+1);
	globalNodePtr[ID+1].createLeaf(start+ 4,4.0);	  
	getLeafBoundsAndComputeAccelerationData(start+8,4,ID+2);
	globalNodePtr[ID+2].createLeaf(start+ 8,4.0);	  
	handleLargeLeaves(ID+3,start+12,items-12);
      }  
  }

  // ==================================================================================================================
  // ==================================================================================================================
  // ==================================================================================================================

  void ParallelQBVHBuilderMortonCode::thread_refit_qbvh(const unsigned int threadID)
  {
    while(1)
      {
	const unsigned int ID = ParallelQBVHBuilder::atomicID.add(1);
	if (ID >= subTrees) break;
	refit_qbvh(subTreeIDs[ID]);
      }
  }

  // ==================================================================================================================
  // ==================================================================================================================
  // ==================================================================================================================

  void ParallelQBVHBuilderMortonCode::refit_qbvh(const unsigned int index)
  {    
    BVHNode &entry = globalNodePtr[index];

    if (unlikely(entry.isLeaf()))    
      {
      }
    else
      {
	const unsigned int children = entry.firstChildID();
	BVHNode *next = &globalNodePtr[children+0];
	prefetch<PFHINT_L1EX>(next + 0);
	prefetch<PFHINT_L1EX>(next + 2);
	prefetch<PFHINT_L2EX>(next + 4);
	prefetch<PFHINT_L2EX>(next + 6);
	const unsigned int items = entry.items();


	MIC_ALIGN Centroid_Scene_AABB childBounds;
	childBounds.reset();

	const AABBExtData e0 = entry.ext_min;
	const AABBExtData e1 = entry.ext_max;
	for (unsigned int i=0;i<items;i++) 
	  {
	    const unsigned int childIndex = children + i;	    	    
	    if (!next[i].isLeaf())
	      refit_qbvh(childIndex);
	    childBounds.extend_scene(upconv4f((float*)&next[i].m_min),
				     upconv4f((float*)&next[i].m_max));
	  }      


	childBounds.storeSceneAABB((float*)&entry);
	entry.ext_min = e0;
	entry.ext_max = e1;
      }    
  }

  // ==================================================================================================================
  // ==================================================================================================================
  // ==================================================================================================================

  void ParallelQBVHBuilderMortonCode::refit_subtrees(const unsigned int index, 
						     const unsigned int threshold, 
						     const unsigned int mode)
  {
    BVHNode &entry = globalNodePtr[index];
    const AABBExtData e0 = entry.ext_min;
    const AABBExtData e1 = entry.ext_max;

    if (!entry.isLeaf())
      {
	const unsigned int children = entry.firstChildID();

	if (e1.t <= threshold)
	  {
	    if (mode == EXTRACT_SUBTREES)
	      {
		if (unlikely(subTrees >= MAX_REFIT_SUBTREES)) FATAL("too many subtrees");
		subTreeIDs[subTrees++] = index;	      	      
		assert(subTrees < MAX_REFIT_SUBTREES);
	      }
	    return;
	  }
      
	if (mode == EXTRACT_SUBTREES)
	  {
	    for (unsigned int i=0;i<entry.items();i++) 
	      {
		const unsigned int childIndex = children + i;	    
		refit_subtrees(childIndex,threshold,mode);	  
	      }
	  }
	else if (mode == REFIT_FROM_SUBTREES)
	  {
	    MIC_ALIGN Centroid_Scene_AABB childBounds;
	    childBounds.reset();

	    for (unsigned int i=0;i<entry.items();i++) 
	      {
		const unsigned int childIndex = children + i;	    
		refit_subtrees(childIndex,threshold,mode);
		childBounds.extend_scene(upconv4f((float*)&globalNodePtr[childIndex].m_min),
					 upconv4f((float*)&globalNodePtr[childIndex].m_max));
	      }      
	    childBounds.storeSceneAABB((float*)&entry);
	    entry.ext_min = e0;
	    entry.ext_max = e1;
	  }
      }
  }

  // ==================================================================================================================
  // ==================================================================================================================
  // ==================================================================================================================


  void ParallelQBVHBuilderMortonCode::build()
  {					
    TIMER(Timer timer);
    TIMER(unsigned long cycles = 0);

    if (unlikely(!triangles)) {
      // should actually never happen
      const mic_f init_node = upconv16f((float*)ParallelQBVHBuilder::initQBVHNode);
      store16f_ngo(&globalNodePtr[0],init_node);
      store16f_ngo(&globalNodePtr[2],init_node);
      store16f_ngo(&globalNodePtr[4],init_node);
      store16f_ngo(&globalNodePtr[6],init_node);
      return;
    }


    numMortonIDBlocks = (triangles+7)/8;
    mortonID[0] = ((MortonID32Bit*)aabb) + 0; 
    mortonID[1] = ((MortonID32Bit*)globalAccelPtr);
    radixTree = (RadixNode*)(((MortonID32Bit*)aabb) + triangles + (16- triangles%16));

    enableTaskStealing = true;

    const unsigned int blocksPerThread = numMortonIDBlocks / NUM_TOTAL_THREADS;


    // ============================================
    // ============================================
    // ============================================
    TIMER(timer.start());  
    Centroid_Scene_AABB &cs_AABB = global_cs_AABB;
    cs_AABB.reset();

    global_cs_AABB.reset();
    splitData.worker_threads = NUM_TOTAL_THREADS;
    splitData.work_items_per_thread = vertices / NUM_TOTAL_THREADS;
    QBVHTaskScheduler::dispatchTask( thread_getSceneBounds );

    TIMER(cycles = (unsigned long)timer.stop());
    TIMER(std::cout << "get scene bounds in " << cycles/1000000.0f << " mcycles " << std::endl);


    // ============================================
    // ============================================
    // ============================================

    // === triangle input array is padded to 16 entry granularity and filled with zeroes !!! ===

    TIMER(timer.start());  
    splitData.work_items_per_thread = blocksPerThread;
    QBVHTaskScheduler::dispatchTask( thread_computeMortonIDs );

    // --- pad to 8 entry mortonID blocks ---
    for (unsigned long i=triangles;i<triangles + (8- triangles%8);i++)
      {
	mortonID[0][i].code  = 0xffffffff;
	mortonID[0][i].index = 0;
      }

    TIMER(cycles = (unsigned long)timer.stop());
    TIMER(std::cout << "create morton codes in " << cycles/1000000.0f << " mcycles " << std::endl << std::flush);

    // ============================================
    // ============================================
    // ============================================

    TIMER(timer.start());  

    splitData.work_items_per_thread = blocksPerThread;
    QBVHTaskScheduler::dispatchTask( thread_radixsort_blocked );
  
    TIMER(cycles = (unsigned long)timer.stop());
    TIMER(const float mtris_per_sec = 1000000000.0f / (float)cycles * (float)triangles / 1000000.0f);
    TIMER(std::cout << "sort morton codes in " << cycles/1000000.0f << " mcycles -> " << mtris_per_sec << " mortonIDs per sec " << std::endl);
  

    // ============================================
    // ============================================
    // ============================================

#if defined(DEBUG)
    for (unsigned long i=1;i<triangles;i++)
      {
	assert(mortonID[0][i].index != 0xffffffff);

	if (mortonID[0][i] < mortonID[0][i-1])
	  {
	    std::cout << "ERROR" << std::endl;
	    DBG_PRINT(i);
	    DBG_PRINT( mortonID[0][i] );
	    DBG_PRINT( mortonID[0][i-1] );
	    FATAL("HERE");
	  }
      }
#endif

    // ============================================
    // ============================================
    // ============================================

    TIMER(timer.start());  
    change01[0].index = 0;
    change01[0].code  = mortonID[0][0].code;
    ParallelQBVHBuilder::atomicID.reset(1);

    splitData.work_items_per_thread = blocksPerThread;
    QBVHTaskScheduler::dispatchTask( thread_getTaskIntervals );

    changes01 = ParallelQBVHBuilder::atomicID.val();

    quicksort_ascending<MortonID32Bit>(change01,1,changes01-1);

    change01[changes01].index = triangles;
    change01[changes01].code = mortonID[0][triangles-1].code;
    changes01++;

    TIMER(DBG_PRINT(changes01));
    assert(changes01 < MAX_TOP_LEVEL_BINS);
    TIMER(cycles = (unsigned long)timer.stop());
    TIMER(std::cout << "get task intervals in " << cycles/1000000.0f << " mcycles" << std::endl);

    // ============================================
    // ============================================
    // ============================================

    TIMER(timer.start());  
    ParallelQBVHBuilder::atomicID.reset(4);

    unsigned long i_index = 0;
    createTopLevel(0,changes01-1,0,change01,i_index);

    TIMER(cycles = (unsigned long)timer.stop());
    TIMER(std::cout << "create top level tree " << cycles/1000000.0f << " mcycles " << std::endl);

    // ============================================
    // ============================================
    // ============================================


    TIMER(timer.start());  

    numActiveLocalWorkQueues =  min(max((changes01+NUM_LOCAL_ATOMIC_BUILD_RECORDS-1) / NUM_LOCAL_ATOMIC_BUILD_RECORDS,(unsigned int)NUM_TOTAL_CORES),(unsigned int)MAX_MIC_CORES);

    if (unlikely(numActiveLocalWorkQueues != NUM_TOTAL_CORES))
      {
	createBuildRecordsFromTaskIntervals();  
      }
    else
      QBVHTaskScheduler::dispatchTask( thread_createBuildRecordsFromTaskIntervals );

    TIMER(cycles = (unsigned long)timer.stop());
    TIMER(std::cout << "transfer build records " << cycles/1000000.0f << " mcycles " << std::endl);

    // ============================================
    // ============================================
    // ============================================

    recursePtr = recurseMortonID;

    TIMER(timer.start());  
    QBVHTaskScheduler::dispatchTask( ParallelQBVHBuilder::thread_build_local );
    TIMER(cycles = (unsigned long)timer.stop());
    TIMER(std::cout << "thread_build " << cycles/1000000.0f << " mcycles " << std::endl);

    // ============================================
    // ============================================
    // ============================================
  
    TIMER(timer.start());  

    qbvh_nodes = atomicID.val() >> 2;

#if 0
    refit_qbvh(0); // === scalar code path ===
#else
    ParallelQBVHBuilder::atomicID.reset(0);
    const unsigned int refit_threshold = (triangles+REFIT_TRIANGLES_PER_THREAD-1) / REFIT_TRIANGLES_PER_THREAD;
    subTrees = 0;
    refit_subtrees(0,refit_threshold,EXTRACT_SUBTREES);

    TIMER(std::cout << "extract subtrees " << (unsigned long)timer.stop()/1000000.0f << " mcycles " << std::endl);

    // ============================================
    // ============================================
    // ============================================
    QBVHTaskScheduler::dispatchTask( ParallelQBVHBuilderMortonCode::thread_refit_qbvh );

    refit_subtrees(0,refit_threshold,REFIT_FROM_SUBTREES);

#endif

    TIMER(cycles = (unsigned long)timer.stop());
    TIMER(std::cout << "refit " << cycles/1000000.0f << " mcycles " << std::endl);

    // ============================================
    // ============================================
    // ============================================
  
  }

  // ==============================================================================================================
  // ==============================================================================================================
  // ==============================================================================================================




  void ParallelQBVHBuilder::check_tree(const bool checkPrimBounds)
  {  
    MIC_ALIGN struct {
      unsigned int node;
    } stackNode[QBVH_MAX_STACK_DEPTH];

    unsigned int nodes = 0;
    unsigned int binaryNodes = 0;
    unsigned int triNodes = 0;
    unsigned int quadNodes = 0;
    unsigned int leaves = 0;
    unsigned int items = 0;
    float trav = 0.0f;
    float isec = 0.0f;
  


    unsigned int equal_parent_child = 0;

    unsigned int sindex = 0;

    for (int i=0;i<1;i++)
      {
	stackNode[i].node = i;
	sindex++;
      }

    unsigned int itemDistribution[32];
    for (unsigned int i=0;i<32;i++)
      itemDistribution[i] = 0;

    DBG_CHECK(DBG_PRINT(globalNodePtr[0]));
    DBG_PRINT(triangles);
    unsigned int *primID = new unsigned int[triangles];
    for (unsigned int i=0;i<triangles;i++) primID[i] = 0;

    nodes += 4;

    while (1) {
      if (unlikely(sindex == 0)) break;
    
      sindex--;
      unsigned int index = stackNode[sindex].node;

      DBG_CHECK(std::cout << "stack pop " << index << std::endl);

      BVHNode &entry = globalNodePtr[index];
      trav += entry.area();

      if (entry.isLeaf())
	{
	  DBG_CHECK(std::cout << "LEAF: offset " << entry.itemListOfs() << " items:" << entry.items() << std::endl);
	  items += entry.items();
	  if(entry.items()<32)
	    itemDistribution[entry.items()]++;
	  else
	    std::cout << "LEAF " << index << " items(" << entry.items() << ") >= 32" << std::endl;

	  isec += entry.area() * entry.items();
	  leaves++;
	  const unsigned int offset = entry.itemListOfs();

	  if(entry.items() > 4) std::cout << "WARNING: " << entry.items() << " items in leaf" << std::endl << std::flush;

	  MIC_ALIGN AABB primitiveBounds;
	  primitiveBounds.setEmpty();
	

	  for (unsigned int i=0;i<entry.items();i++)
	    {	    
	      if (checkPrimBounds) primitiveBounds.extend(aabb[offset+i]);
	      primID[offset+i]++;
	    }
#if 1
	  if (checkPrimBounds)
	    if(!entry.enclose(primitiveBounds))
	      {	    
		DBG_PRINT(entry);
		DBG_PRINT(primitiveBounds);
		std::cout << "offset " << offset << " items " << entry.items() << std::endl;
		//exit(0);
	      }
#endif
	}
      else
	{
	  DBG_CHECK(std::cout << "NODE: offset " << entry.firstChildID() << std::endl);

	  assert(sindex + 4 < QBVH_MAX_STACK_DEPTH);
	  const unsigned int children = entry.firstChildID();
	  unsigned int n = 0;

	  DBG_CHECK(DBG_PRINT(entry.items()));
	  nodes += entry.items(); 
	  //nodes += 4; 

	  if (entry.items() == 2) 
	    binaryNodes++;
	  else if (entry.items() == 3) 
	    triNodes++;
	  else if (entry.items() == 4) 
	    quadNodes++;
	  else
	    {
	      DBG_PRINT(entry.items());
	      std::cout << index << " number of nodes " << entry.items() << std::endl; 
	    } 

	  unsigned int child_node_leaves = 0;
	  for (unsigned int i=0;i<entry.items();i++) 
	    {
	      n++;
	      const unsigned int childIndex = children + i;
	    
	      if (globalNodePtr[childIndex].isLeaf())
		child_node_leaves++;

	      if(!entry.enclose(globalNodePtr[childIndex]))
		{
		  std::cout << "Enclosing error" << std::endl;
		  std::cout << "parent " << index << " -> " << entry << std::endl;
		  std::cout << "child  " << childIndex << " -> " << globalNodePtr[childIndex] << std::endl;
		  exit(0);
		}

	      if (equal(entry,globalNodePtr[childIndex]))
		equal_parent_child++;

	      stackNode[sindex+i].node = childIndex;
	    }      

	  sindex+=n;
	  DBG_CHECK(DBG_PRINT(sindex));
	}
    
    }

    for (unsigned int i=0;i<triangles;i++) 
      if (primID[i] == 0 || primID[i] > 1)
	{
	  std::cout << "primitive error at " << i << " " << primID[i] << std::endl;
	  exit(0);
	}

    int sum = 0;
  
    for (unsigned int i=0;i<32;i++)
      sum += itemDistribution[i];

    int part = 0;
    for (unsigned int i=0;i<8;i++)
      {
	part += itemDistribution[i];
	std::cout << i << " [" << 100.0f * itemDistribution[i] / sum << ", " << 100.0f * part / sum << " p-sum " << part << "] ";
      }
    std::cout << std::endl;

    part = 0;
    for (unsigned int i=0;i<=4;i++)
      {
	part += i*itemDistribution[i];
      }

    std::cout << "leaf/item util        " << 100.0f * (float)part / (sum*4) << std::endl;
    std::cout << "equal_parent_child    " << equal_parent_child << std::endl;

    DBG_PRINT(globalNodePtr[0].area());
    DBG_PRINT(isec);
    DBG_PRINT(trav);
    DBG_PRINT(isec / globalNodePtr[0].area());
    DBG_PRINT(trav / globalNodePtr[0].area());

    
    DBG_PRINT(nodes);
    std::cout << "binaryNodes " << binaryNodes << " [" << binaryNodes*2 << "]" << std::endl;
    std::cout << "triNodes    " << triNodes << " [" << triNodes*4 << "]" << std::endl;
    std::cout << "quadNodes   " << quadNodes << " [" << quadNodes*4 << "]" << std::endl;
    std::cout << "BVH size " << (float)sizeof(BVHNode)*nodes / 1024.0f << " KB" << std::endl;
    const float node_util = (float)(2*binaryNodes+3*triNodes+4*quadNodes) / (binaryNodes+triNodes+quadNodes);
    std::cout << "quad node utilization " << node_util <<  " " << 100.0f * node_util / 4 << std::endl;
  
    DBG_PRINT(leaves);
    DBG_PRINT(items);      
    DBG_PRINT(trav);
    DBG_PRINT(isec);

    delete [] primID;
  }


};

