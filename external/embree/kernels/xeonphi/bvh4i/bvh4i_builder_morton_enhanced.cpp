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
#include "bvh4i/bvh4i_builder_morton_enhanced.h"
#include "limits.h"


#define SINGLE_THREADED_BUILD_THRESHOLD  (MAX_MIC_THREADS*8)

#define DIAG_FACTOR 0.15f
#define MAX_REBUILD_NODES numPrimitives
#define PROFILE_ITERATIONS 20

#define L1_PREFETCH_ITEMS 2
#define L2_PREFETCH_ITEMS 16

#define TIMER(x) 

//#define PROFILE
//#define TEST_BUILD_PERFORMANCE

#define DBG(x) 

namespace embree 
{
  AtomicMutex mtx;
  __aligned(64) static double dt = 0.0f;

  BVH4iBuilderMortonEnhanced::BVH4iBuilderMortonEnhanced (BVH4i* _bvh, BuildSource* _source, void* _geometry)
    : BVH4iBuilderMorton(_bvh,_source,_geometry){}

  void BVH4iBuilderMortonEnhanced::build(size_t threadIndex, size_t threadCount) 
  {
    DBG(PING);

    if (g_verbose >= 2) {
      std::cout << "building BVH4i with Enhanced Morton builder (MIC)... " << std::endl << std::flush;
    }

    initEncodingAllocateData(threadCount);
    LockStepTaskScheduler::init(TaskScheduler::getNumThreads()); 

    DBG_PRINT( numPrimitives );
#if defined(PROFILE)
    std::cout << "STARTING PROFILE MODE" << std::endl << std::flush;

    double dt_min = pos_inf;
    double dt_avg = 0.0f;
    double dt_max = neg_inf;
    size_t iterations = PROFILE_ITERATIONS;
    for (size_t i=0; i<iterations; i++) 
      {
	TaskScheduler::executeTask(threadIndex,threadCount,_build_parallel_morton_enhanced,this,TaskScheduler::getNumThreads(),"build_parallel_morton_enhanced");

	dt_min = min(dt_min,dt);
	dt_avg = dt_avg + dt;
	dt_max = max(dt_max,dt);
      }
    dt_avg /= double(iterations);

    std::cout << "[DONE]" << std::endl;
    std::cout << "  min = " << 1000.0f*dt_min << "ms (" << numPrimitives/dt_min*1E-6 << " Mtris/s)" << std::endl;
    std::cout << "  avg = " << 1000.0f*dt_avg << "ms (" << numPrimitives/dt_avg*1E-6 << " Mtris/s)" << std::endl;
    std::cout << "  max = " << 1000.0f*dt_max << "ms (" << numPrimitives/dt_max*1E-6 << " Mtris/s)" << std::endl;
    std::cout << BVH4iStatistics(bvh).str();

#else
    DBG(DBG_PRINT(numPrimitives));


    if (likely(numPrimitives > SINGLE_THREADED_BUILD_THRESHOLD && TaskScheduler::getNumThreads() > 1))
      {
	DBG(std::cout << "PARALLEL BUILD" << std::endl << std::flush);
	TaskScheduler::executeTask(threadIndex,threadCount,_build_parallel_morton_enhanced,this,TaskScheduler::getNumThreads(),"build_parallel");
      }
    else
      {
	/* number of primitives is small, just use single threaded mode */
	if (likely(numPrimitives > 0))
	  {
	    DBG(std::cout << "SERIAL BUILD" << std::endl << std::flush);
	    build_parallel_morton_enhanced(0,1,0,0,NULL);
	  }
	else
	  {
	    DBG(std::cout << "EMPTY SCENE BUILD" << std::endl << std::flush);
	    /* handle empty scene */
	    for (size_t i=0;i<4;i++)
	      bvh->qbvh[0].setInvalid(i);
	    for (size_t i=0;i<4;i++)
	      bvh->qbvh[1].setInvalid(i);
	    bvh->qbvh[0].lower[0].child = BVH4i::NodeRef(128);
	    bvh->root = bvh->qbvh[0].lower[0].child; 
	    bvh->bounds = BBox3f(*(Vec3fa*)&bvh->qbvh->lower[0],*(Vec3fa*)&bvh->qbvh->upper[0]);	    
	  }
      }

    if (g_verbose >= 2) {
      double perf = numPrimitives/dt*1E-6;
      std::cout << "[DONE] " << 1000.0f*dt << "ms (" << perf << " Mtris/s), primitives " << numPrimitives << std::endl;
      std::cout << BVH4iStatistics(bvh).str();
    }
#endif

  }

  bool splitSAH(PrimRef * __restrict__ const prims, BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild)
  {
    const unsigned int items = current.items();

    /* mark as leaf if leaf threshold reached */
    if (items <= BVH4i::N) {
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

    if (leftChild.items()  <= BVH4i::N) leftChild.createLeaf();
    if (rightChild.items() <= BVH4i::N) rightChild.createLeaf();	
    return true;
  }



  void BVH4iBuilderMortonEnhanced::extractTopLevelTree(const size_t index,
						       const Vec3fa &root_diag,
						       BVHNode *__restrict__ const local_node,
						       size_t &nodes)
  {
    BVHNode &entry = this->node[index];

    if (!entry.isLeaf())
      {
	const unsigned int children = entry.firstChildID();

	const Vec3fa diag = entry.upper - entry.lower;
	if (diag.x < DIAG_FACTOR * root_diag.x &&
	    diag.y < DIAG_FACTOR * root_diag.y &&
	    diag.z < DIAG_FACTOR * root_diag.z)
	  {
	    if (unlikely(nodes >= MAX_REBUILD_NODES)) FATAL("too many subtrees");
	    local_node[nodes++] = entry;
	    assert(nodes < MAX_REBUILD_NODES);
	    return;
	  }
	else
	  {
	    for (unsigned int i=0;i<entry.items();i++) 
	      {
		const unsigned int childIndex = children + i;	    
		extractTopLevelTree(childIndex,root_diag,local_node,nodes);	  
	      }
	  }
      }
    else
      {
	if (unlikely(nodes >= MAX_REBUILD_NODES)) FATAL("too many subtrees");
	local_node[nodes++] = entry;
	assert(nodes < MAX_REBUILD_NODES);
      }     
  }


  void BVH4iBuilderMortonEnhanced::buildTopLevelSAHTree(BVHNode &parent,
							BuildRecord &current,
							BVHNode *__restrict__ local_node)
  {
    DBG(std::cout << std::endl; PING);
    DBG(DBG_PRINT(parent));
    DBG(DBG_PRINT(current));
#ifdef DEBUG
    {
      PrimRef tmp;
      tmp.lower = current.bounds.geometry.lower;
      tmp.upper = current.bounds.geometry.upper;
      for (size_t i=0;i<current.items();i++)
	assert(subset(*(PrimRef*)&local_node[current.begin+i],tmp) == true);
    }
#endif      

    BuildRecord record[BVH4i::N];
    BuildRecord left, right;

    unsigned int numSplits = 1;
    record[0] = current;

    while(1)
      {
	bool couldSplit = false;
	if (numSplits < BVH4i::N)
	  {
	    int index = -1;
	    float maxArea = -1.0f;
	    for (unsigned int i=0;i<numSplits;i++)
	      {
		if (!record[i].isLeaf())
		  {
		    assert(record[i].sceneArea() > 0.0f);

		    if (record[i].sceneArea() >= maxArea)
		      {
			maxArea = record[i].sceneArea();
			index = i;
		      }
		  }
	      }

	    if (index != -1)
	      {
		assert(index < BVH4i::N);
		bool s = splitSAH((PrimRef*)local_node,record[index],left,right);  

		assert(numSplits > 0);

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
		    DBG(std::cout << "became leaf" << std::endl);
		    continue;
		  }
	      }
	  }
	// ==================================================
	if (!couldSplit) break;
      }

    DBG(std::cout << "DETERMINED ALL SPLITS" << std::endl);


    DBG(DBG_PRINT(numSplits));

    // could not generate any split, propagate leaf back to parent
    if ( numSplits == 1 )
      {
	current = record[0];
	DBG(std::cout << "CREATING LOCAL NODE LEAF WITH " << current.items() << " NODES" << std::endl);
	assert(current.isLeaf());
	if(current.items() > BVH4i::N) 
	  {
	    std::cout << "WARNING UNSPLITABLE LEAF " << std::endl;
	    FATAL("SHOULD NOT HAPPEN FOR TOP-LEVEL TREE");
	    DBG_PRINT(current.items());
	    record[1] = record[0];
	    record[0].end   = record[0].begin + 4;
	    record[1].begin = record[0].begin + 4;
	    numSplits = 2;
	  }
	else
	  {
	    if (unlikely(current.items() == 1))
	      {
		parent = local_node[current.begin];
		return;
	      }
	    assert(current.items() > 1);

	    const unsigned int currentIndex = this->atomicID.add(BVH4i::N);
	    if (unlikely(currentIndex >= this->numAllocatedNodes))
	      {
		DBG_PRINT(this->numAllocatedNodes);
		FATAL("not enough nodes allocated");
	      }
      
	    /* init used/unused nodes */
	    const mic_f init_node = load16f((float*)BVH4i::initQBVHNode);
	    store16f((float*)&node[currentIndex+0],init_node);
	    store16f((float*)&node[currentIndex+2],init_node);

	    DBG(DBG_PRINT(currentIndex));

	    parent.createNode(currentIndex,current.items());
	    //current.childrenID = currentIndex;

	    for (size_t i=0;i<current.items();i++)
	      {
		assert(subset(*(PrimRef*)&local_node[current.begin+i],*(PrimRef*)&parent));
		this->node[currentIndex+i] = local_node[current.begin+i];
	      }

	    DBG(DBG_PRINT(parent));
	    return;
	  }
      }

    assert(numSplits >= 2 && numSplits <= BVH4i::N);

    // ==== aquire next four nodes ====
    const unsigned int currentIndex = this->atomicID.add(BVH4i::N);
    if (unlikely(currentIndex >= this->numAllocatedNodes))
      {
	DBG_PRINT(this->numAllocatedNodes);
	FATAL("not enough nodes allocated");
      }

    const mic_f init_node = load16f((float*)BVH4i::initQBVHNode);
    store16f((float*)&node[currentIndex+0],init_node);
    store16f((float*)&node[currentIndex+2],init_node);

    parent.createNode(currentIndex,numSplits);
    for (unsigned int i=0;i<numSplits;i++)
      {
	DBG(DBG_PRINT(currentIndex+i));
        this->node[currentIndex+i].lower = record[i].bounds.geometry.lower;
        this->node[currentIndex+i].upper = record[i].bounds.geometry.upper;
	buildTopLevelSAHTree(this->node[currentIndex+i],record[i],local_node);
      }

  }


  void BVH4iBuilderMortonEnhanced::build_parallel_morton_enhanced(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event) 
  {
    DBG(PING);
    /* initialize thread state */
    initThreadState(threadIndex,threadCount);
    
    /* let all thread except for control thread wait for work */
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

    // ===================          
    // === start timer ===
    // ===================
    
    TIMER(double msec);
    TIMER(msec = getSeconds());

    build_main(threadIndex,threadCount);

    bvh->accel = this->accel;
    bvh->qbvh  = (BVH4i::Node*)this->node;

    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "build morton tree " << 1000. * msec << " ms" << std::endl << std::flush);

    // ==============================        
    // === extract top level tree ===
    // ==============================


    TIMER(msec = getSeconds());

    const Vec3fa rootDiag = this->node[0].upper - this->node[0].lower;
    BVHNode * __restrict__ local_node = (BVHNode*)morton;
    size_t nodes = 0;
    extractTopLevelTree(0,rootDiag,local_node,nodes);

    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "extract top-level nodes " << 1000. * msec << " ms" << std::endl << std::flush);

    TIMER(DBG_PRINT(nodes));

    // ==============================        
    // === rebuild top level tree ===
    // ==============================

    TIMER(msec = getSeconds());    
    BuildRecord topLevelBuildRecord;
    topLevelBuildRecord.init(this->global_bounds,0,nodes);
    buildTopLevelSAHTree(this->node[0],topLevelBuildRecord,local_node);

    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "rebuild top-level " << 1000. * msec << " ms" << std::endl << std::flush);
	    
    // ===================================        
    // === convert to optimized layout ===
    // ===================================

    this->numNodes = this->atomicID >> 2;

    TIMER(msec = getSeconds());    

    LockStepTaskScheduler::dispatchTask( task_convertToSOALayout, this, threadIndex, threadCount );

    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "task_convertToSOALayout " << 1000. * msec << " ms" << std::endl << std::flush);


    bvh->root = bvh->qbvh[0].lower[0].child; 
    bvh->bounds = BBox3f(*(Vec3fa*)&bvh->qbvh->lower[0],*(Vec3fa*)&bvh->qbvh->upper[0]);
    
    // ==================          
    // === stop timer ===
    // ==================

    if (g_verbose >= 2) {
      double t1 = getSeconds();
      double perf = numPrimitives/(t1-t0)*1E-6;
      std::cout << "[DONE] " << t1-t0 << "sec (" << perf << " Mtris/s)" << std::endl;
      //
    }

    LockStepTaskScheduler::releaseThreads(threadCount);

    /* stop measurement */
#if !defined(PROFILE)
    if (g_verbose >= 2) 
#endif
      dt = getSeconds()-t0;

  };

};
