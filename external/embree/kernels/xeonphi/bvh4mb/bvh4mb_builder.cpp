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

#include "bvh4mb/bvh4mb_builder.h"

namespace embree
{
#define DBG(x) 
#define TIMER(x) 

#define L1_PREFETCH_ITEMS 2
#define L2_PREFETCH_ITEMS 16

#define GENERATE_SUBTREES_MAX_TREE_DEPTH 6
#define SERIAL_REFIT_THRESHOLD 1024

  Builder* BVH4mbBuilder::create (void* accel, BuildSource* source, void* geometry, size_t mode ) 
  { 
    Builder* builder = new BVH4mbBuilder((BVH4mb*)accel,source,geometry);
    return builder;
  }

  void BVH4mbBuilder::printBuilderName()
  {
    std::cout << "building BVH4mb with binned SAH builder (MIC) ... " << std::endl;    
  }


  size_t BVH4mbBuilder::getNumPrimitives()
  {
    if (scene->numTriangleMeshes2 == 0) return 0;

    /* count total number of triangles */
    size_t primitives = 0;       
    for (size_t i=0;i<scene->size();i++)
      {
	if (unlikely(scene->get(i) == NULL)) continue;
	if (unlikely((scene->get(i)->type != TRIANGLE_MESH))) continue;
	if (unlikely(!scene->get(i)->isEnabled())) continue;

	const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(i);
	if (unlikely(mesh->numTimeSteps == 1)) continue;

	primitives += mesh->numTriangles;
      }
    return primitives;	  
  }

  void BVH4mbBuilder::computePrimRefs(const size_t threadIndex, const size_t threadCount)
  {
    LockStepTaskScheduler::dispatchTask( task_computePrimRefsTrianglesMB, this, threadIndex, threadCount );
  }


  void BVH4mbBuilder::computePrimRefsTrianglesMB(const size_t threadID, const size_t numThreads) 
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
      if (unlikely(mesh->numTimeSteps == 1)) continue;

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
      if (unlikely(mesh->numTimeSteps == 1)) continue;

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


  void BVH4mbBuilder::allocateData(const size_t threadCount, const size_t totalNumPrimitives)
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
	allocateMemoryPools(numPrims,numNodes,sizeof(BVH4i::Node),sizeof(BVH4mb::Triangle01));
      }
  }

  __forceinline void computeAccelerationDataMB(const unsigned int &geomID,
					     const unsigned int &primID,     
					     const Scene *__restrict__ const scene,
					     BVH4mb::Triangle01 * __restrict__ const acc)
  {
    const TriangleMeshScene::TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
    const TriangleMeshScene::TriangleMesh::Triangle & tri = mesh->triangle(primID);

    const mic_i pID(primID);
    const mic_i gID(geomID);

    const float *__restrict__ const vptr0_t0 = (float*)&mesh->vertex(tri.v[0]);
    const float *__restrict__ const vptr1_t0 = (float*)&mesh->vertex(tri.v[1]);
    const float *__restrict__ const vptr2_t0 = (float*)&mesh->vertex(tri.v[2]);

    prefetch<PFHINT_L1>(vptr1_t0);
    prefetch<PFHINT_L1>(vptr2_t0);

    const mic_f v0_t0 = broadcast4to16f(vptr0_t0); 
    const mic_f v1_t0 = broadcast4to16f(vptr1_t0);
    const mic_f v2_t0 = broadcast4to16f(vptr2_t0);

    const mic_f tri_accel_t0 = initTriangle1(v0_t0,v1_t0,v2_t0,gID,pID,mic_i(mesh->mask));

    store16f_ngo(&acc->t0,tri_accel_t0);

    if ((int)mesh->numTimeSteps == 1)
      {
	store16f_ngo(&acc->t1,tri_accel_t0);
      }
    else
      {
	assert( (int)mesh->numTimeSteps == 2 );
	const float *__restrict__ const vptr0_t1 = (float*)&mesh->vertex(tri.v[0],1);
	const float *__restrict__ const vptr1_t1 = (float*)&mesh->vertex(tri.v[1],1);
	const float *__restrict__ const vptr2_t1 = (float*)&mesh->vertex(tri.v[2],1);
	
	const mic_f v0_t1 = broadcast4to16f(vptr0_t1); 
	const mic_f v1_t1 = broadcast4to16f(vptr1_t1);
	const mic_f v2_t1 = broadcast4to16f(vptr2_t1);

	const mic_f tri_accel_t1 = initTriangle1(v0_t1,v1_t1,v2_t1,gID,pID,mic_i(mesh->mask));

	store16f_ngo(&acc->t1,tri_accel_t1);
      }
  }

  void BVH4mbBuilder::createTriangle01AccelMB(const size_t threadID, const size_t numThreads)
  {
    const size_t startID = (threadID+0)*numPrimitives/numThreads;
    const size_t endID   = (threadID+1)*numPrimitives/numThreads;

    BVH4mb::Triangle01 * __restrict__  acc  = (BVH4mb::Triangle01 *)accel + startID;
    const PrimRef* __restrict__  bptr = prims + startID;

    for (size_t j=startID; j<endID; j++, bptr++, acc++)
      {
	prefetch<PFHINT_NT>(bptr + L1_PREFETCH_ITEMS);
	prefetch<PFHINT_L2>(bptr + L2_PREFETCH_ITEMS);
	assert(bptr->geomID() < scene->size() );
	assert(bptr->primID() < scene->get( bptr->geomID() )->numPrimitives );

	computeAccelerationDataMB(bptr->geomID(),bptr->primID(),scene,acc);
      }
  }

  void BVH4mbBuilder::createAccel(const size_t threadIndex, const size_t threadCount)
  {
    DBG(PING);
    LockStepTaskScheduler::dispatchTask( task_createTriangle01AccelMB, this, threadIndex, threadCount );   
  }

  void BVH4mbBuilder::refit(const size_t index)
  {   
    BVHNode& entry = node[index];
    if (unlikely(entry.isLeaf()))
      {
	unsigned int accel_entries = entry.items();
	unsigned int accel_offset  = entry.itemListOfs();
	BBox3f leaf_bounds = empty;
	BVH4mb::Triangle01* accelMB = (BVH4mb::Triangle01*)accel + accel_offset;
	for (size_t i=0;i<accel_entries;i++)
	  {
	    leaf_bounds.extend( accelMB[i].t1.bounds() );
	  }

	*(BBox3f*)&node[index+4] = leaf_bounds;
	return;
      }

    const size_t childrenID = entry.firstChildID();
    const size_t items    = entry.items();
    BBox3f* next = (BBox3f*)&node[childrenID+4];
    
    /* init second node */
    const mic_f init_node = load16f((float*)BVH4i::initQBVHNode);
    store16f_ngo(next + 0,init_node);
    store16f_ngo(next + 2,init_node);

    BBox3f parentBounds = empty;
    for (size_t i=0; i<items; i++) 
    {
      const size_t childIndex = childrenID + i;	    	    
      refit(childIndex);
    }      


    for (size_t i=0; i<items; i++) 
      parentBounds.extend( next[i] );

    *(BBox3f*)&node[index+4] = parentBounds;    

  }    

  void BVH4mbBuilder::check_tree(const unsigned index)
  {
    BVHNode& entry = node[index];
    if (unlikely(entry.isLeaf()))
      {
	unsigned int accel_entries = entry.items();
	unsigned int accel_offset  = entry.itemListOfs();
	BBox3f leaf_bounds = empty;
	BVH4mb::Triangle01* accelMB = (BVH4mb::Triangle01*)accel + accel_offset;
	for (size_t i=0;i<accel_entries;i++)
	  leaf_bounds.extend( accelMB[i].t1.bounds() );
	if (leaf_bounds != *(BBox3f*)&node[index+4])
	  {
	    DBG_PRINT(leaf_bounds);
	    DBG_PRINT(*(BBox3f*)&node[index+4]);
	    DBG_PRINT(index);
	    FATAL("LEAF");
	  }
      }
    else
      {
	const size_t childrenID = entry.firstChildID();
	const size_t items    = entry.items();
	BBox3f* next = (BBox3f*)&node[childrenID+4];

	BBox3f bounds = empty;
	for (size_t i=0; i<items; i++) 
	  bounds.extend ( next[i]  );

	if (index != 0)
	  if ( bounds != *(BBox3f*)&node[index+4])
	    {
	      DBG_PRINT(bounds);
	      DBG_PRINT(*(BBox3f*)&node[index+4]);
	      DBG_PRINT(index);
	    FATAL("NODE");
	    }

	for (size_t i=0; i<items; i++) 
	  {
	    const size_t childIndex = childrenID + i;	    	    
	    check_tree(childIndex);
	  }
      }
    
  }

  BBox3f BVH4mbBuilder::refit_subtree(const size_t index)
  {
    BBox3f local[4];

    BVHNode& entry = node[index];
    if (unlikely(entry.isLeaf()))
      {
	unsigned int accel_entries = entry.items();
	unsigned int accel_offset  = entry.itemListOfs();
	BBox3f leaf_bounds = empty;
	BVH4mb::Triangle01* accelMB = (BVH4mb::Triangle01*)accel + accel_offset;

	for (size_t i=0;i<accel_entries;i++)
	  prefetch<PFHINT_L1>(&accelMB[i].t1);

	for (size_t i=0;i<accel_entries;i++)
	  leaf_bounds.extend( accelMB[i].t1.bounds() );
	return leaf_bounds;
      }

    prefetch<PFHINT_L1>(local + 0);
    prefetch<PFHINT_L1>(local + 2);

    const size_t childrenID = entry.firstChildID();
    const size_t items    = entry.items();
    BBox3f* next = (BBox3f*)&node[childrenID+4];


    /* init second node */
    const mic_f init_node = load16f((float*)BVH4i::initQBVHNode);
    store16f(local + 0,init_node);
    store16f(local + 2,init_node);


    BBox3f bounds = empty;
    for (size_t i=0; i<items; i++) 
    {
      const size_t childIndex = childrenID + i;	    	    
      BBox3f childBounds = refit_subtree(childIndex);
      local[i] = childBounds;
      bounds.extend ( childBounds );
    }      

    store16f_ngo(next + 0,load16f(&local[0]));
    store16f_ngo(next + 2,load16f(&local[2]));
    
    return bounds;
    
  }




  void BVH4mbBuilder::generate_subtrees(const size_t index,const size_t depth, size_t &subtrees)
  {
    BVHNode& entry = node[index];

    if (depth == GENERATE_SUBTREES_MAX_TREE_DEPTH || entry.isLeaf())
      {
	unsigned int *subtrees_array = (unsigned int*)prims;
	subtrees_array[subtrees++] = index;
	return;
      }

    const size_t childrenID = entry.firstChildID();
    const size_t items      = entry.items();

    for (size_t i=0; i<items; i++) 
      {
	const size_t childIndex = childrenID + i;	    	    
	generate_subtrees(childIndex,depth+1,subtrees);
      }      
  }


  BBox3f BVH4mbBuilder::refit_toplevel(const size_t index,const size_t depth)
  {
    BVHNode& entry = node[index];

    if (depth == GENERATE_SUBTREES_MAX_TREE_DEPTH || entry.isLeaf())
      {
	return *(BBox3f*)&node[index+4];
      }

    const size_t childrenID = entry.firstChildID();
    const size_t items    = entry.items();
    BBox3f* next = (BBox3f*)&node[childrenID+4];

    __aligned(64) BBox3f local[4];

    /* init second node */
    const mic_f init_node = load16f((float*)BVH4i::initQBVHNode);
    store16f(local + 0,init_node);
    store16f(local + 2,init_node);


    BBox3f bounds = empty;
    for (size_t i=0; i<items; i++) 
    {
      const size_t childIndex = childrenID + i;	    	    
      BBox3f childBounds = refit_toplevel(childIndex,depth+1);
      local[i] = childBounds;
      bounds.extend ( childBounds );
    }      

    store16f_ngo(next + 0,load16f(&local[0]));
    store16f_ngo(next + 2,load16f(&local[2]));
    
    return bounds;
  }

  __forceinline void convertToBVH4MBLayout(BVHNode *__restrict__ const bptr)
  {
    const mic_i box01 = load16i((int*)(bptr + 0));
    const mic_i box23 = load16i((int*)(bptr + 2));

    const mic_i box_min01 = permute<2,0,2,0>(box01);
    const mic_i box_max01 = permute<3,1,3,1>(box01);

    const mic_i box_min23 = permute<2,0,2,0>(box23);
    const mic_i box_max23 = permute<3,1,3,1>(box23);
    const mic_i box_min0123 = select(0x00ff,box_min01,box_min23);
    const mic_i box_max0123 = select(0x00ff,box_max01,box_max23);

    const mic_m min_d_mask = bvhLeaf(box_min0123) != mic_i::zero();
    const mic_i childID    = bvhChildID(box_min0123)>>2;
    const mic_i min_d_node = qbvhCreateNode(childID,mic_i::zero());
    const mic_i min_d_leaf = ((box_min0123 ^ BVH_LEAF_MASK)<<1) | QBVH_LEAF_MASK; // * 2 as accel size is 128 bytes now
    const mic_i min_d      = select(min_d_mask,min_d_leaf,min_d_node);
    const mic_i bvh4_min   = select(0x7777,box_min0123,min_d);
    const mic_i bvh4_max   = box_max0123;
    store16i_nt((int*)(bptr + 0),bvh4_min);
    store16i_nt((int*)(bptr + 2),bvh4_max);
  }

  void BVH4mbBuilder::convertToSOALayoutMB(const size_t threadID, const size_t numThreads)
  {
    const size_t startID = (threadID+0)*numNodes/numThreads;
    const size_t endID   = (threadID+1)*numNodes/numThreads;

    BVHNode  * __restrict__  bptr = ( BVHNode*)node + startID*4;

    BVH4i::Node * __restrict__  qptr = (BVH4i::Node*)node + startID;

    for (unsigned int n=startID;n<endID;n++,qptr++,bptr+=4) 
      {
	prefetch<PFHINT_L1EX>(bptr+4);
	prefetch<PFHINT_L2EX>(bptr+4*4);
	convertToBVH4MBLayout(bptr);
	evictL1(bptr);
      }
  }


  void BVH4mbBuilder::refitBVH4MB(const size_t threadID, const size_t numThreads)
  {
    const size_t startID = (threadID+0)*subtrees/numThreads;
    const size_t endID   = (threadID+1)*subtrees/numThreads;

    unsigned int *subtrees_array = (unsigned int*)prims;

    while(1)
      {
	unsigned int ID = atomicID.inc();
	if (ID >= subtrees) break;
	const unsigned int index = subtrees_array[ID];
	BBox3f bounds = refit_subtree(index);
	*(BBox3f*)&node[index+4] = bounds;		
      }
  }

  void BVH4mbBuilder::convertQBVHLayout(const size_t threadIndex, const size_t threadCount)
  {
    TIMER(double msec = 0.0);


    if (numPrimitives < SERIAL_REFIT_THRESHOLD)
      {
	refit(0);
      }
    else
      {
	TIMER(msec = getSeconds());
	// ------------------------
	atomicID.reset(0);
	subtrees = 0;
	generate_subtrees(0,0,subtrees);
	// ------------------------
	TIMER(msec = getSeconds()-msec);    
	TIMER(std::cout << "generate subtrees " << 1000. * msec << " ms" << std::endl << std::flush);

	DBG(DBG_PRINT(subtrees));

	TIMER(msec = getSeconds());
	// ------------------------
	LockStepTaskScheduler::dispatchTask( task_refitBVH4MB, this, threadIndex, threadCount );    
	// ------------------------
	TIMER(msec = getSeconds()-msec);    
	TIMER(std::cout << "refit subtrees " << 1000. * msec << " ms" << std::endl << std::flush);

	TIMER(msec = getSeconds());
	// ------------------------    
	refit_toplevel(0,0);
	// ------------------------
	TIMER(msec = getSeconds()-msec);    
	TIMER(std::cout << "refit toplevel " << 1000. * msec << " ms" << std::endl << std::flush);

      }

#if defined(DEBUG)
    std::cout << "checking tree..." << std::flush;
    check_tree(0);
    std::cout << "done" << std::endl << std::flush;
#endif


    TIMER(msec = getSeconds());

    LockStepTaskScheduler::dispatchTask( task_convertToSOALayoutMB, this, threadIndex, threadCount );    

    TIMER(msec = getSeconds()-msec);    
    TIMER(std::cout << "convert " << 1000. * msec << " ms" << std::endl << std::flush);

  }

}
