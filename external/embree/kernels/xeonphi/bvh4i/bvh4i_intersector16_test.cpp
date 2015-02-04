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

#include "bvh4i_intersector16_test.h"
#include "bvh4i_leaf_intersector.h"

namespace embree
{
  namespace isa
  {
    
    static unsigned int BVH4I_LEAF_MASK = BVH4i::leaf_mask; // needed due to compiler efficiency bug
    static unsigned int M_LANE_7777 = 0x7777;               // needed due to compiler efficiency bug

    __forceinline void convertSOA4ftoAOS4f(const mic_f &x,
					   const mic_f &y,
					   const mic_f &z,
					   const mic_f &w,
					   Vec3fa *__restrict__ const dest)
    {
      mic_f r0 = mic_f::undefined();
      r0 = uload16f_low(0x1111,&x[0],r0);
      r0 = uload16f_low(0x2222,&y[0],r0);
      r0 = uload16f_low(0x4444,&z[0],r0);
      r0 = uload16f_low(0x8888,&w[0],r0);

      mic_f r1 = mic_f::undefined();
      r1 = uload16f_low(0x1111,&x[4],r1);
      r1 = uload16f_low(0x2222,&y[4],r1);
      r1 = uload16f_low(0x4444,&z[4],r1);
      r1 = uload16f_low(0x8888,&w[4],r1);

      mic_f r2 = mic_f::undefined();
      r2 = uload16f_low(0x1111,&x[8],r2);
      r2 = uload16f_low(0x2222,&y[8],r2);
      r2 = uload16f_low(0x4444,&z[8],r2);
      r2 = uload16f_low(0x8888,&w[8],r2);

      mic_f r3 = mic_f::undefined();
      r3 = uload16f_low(0x1111,&x[12],r3);
      r3 = uload16f_low(0x2222,&y[12],r3);
      r3 = uload16f_low(0x4444,&z[12],r3);
      r3 = uload16f_low(0x8888,&w[12],r3);

      store16f(&dest[ 0],r0);
      store16f(&dest[ 4],r1);
      store16f(&dest[ 8],r2);
      store16f(&dest[12],r3);
    }

    // ============================================================================================
    // ============================================================================================
    // ============================================================================================

    template<typename LeafIntersector, bool ENABLE_COMPRESSED_BVH4I_NODES>
    void BVH4iIntersector16Test<LeafIntersector,ENABLE_COMPRESSED_BVH4I_NODES>::intersect(mic_i* valid_i, BVH4i* bvh, Ray16& ray16)
    {
      /* near and node stack */
      __aligned(64) float   stack_dist[3*BVH4i::maxDepth+1];
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];


      /* setup */
      const mic_m m_valid    = *(mic_i*)valid_i != mic_i(0);
      const mic3f rdir16     = rcp_safe(ray16.dir);
      const mic_f inf        = mic_f(pos_inf);
      const mic_f zero       = mic_f::zero();

      store16f(stack_dist,inf);

      const Node      * __restrict__ nodes = (Node     *)bvh->nodePtr();
      const Triangle1 * __restrict__ accel = (Triangle1*)bvh->triPtr();

      stack_node[0] = BVH4i::invalidNode;
      long rayIndex = -1;
      while((rayIndex = bitscan64(rayIndex,toInt(m_valid))) != BITSCAN_NO_BIT_SET_64)	    
        {
	  stack_node[1] = bvh->root;
	  size_t sindex = 2;

	  const mic_f org_xyz      = loadAOS4to16f(rayIndex,ray16.org.x,ray16.org.y,ray16.org.z);
	  const mic_f dir_xyz      = loadAOS4to16f(rayIndex,ray16.dir.x,ray16.dir.y,ray16.dir.z);
	  const mic_f rdir_xyz     = loadAOS4to16f(rayIndex,rdir16.x,rdir16.y,rdir16.z);
	  const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
	  const mic_f min_dist_xyz = broadcast1to16f(&ray16.tnear[rayIndex]);
	  mic_f       max_dist_xyz = broadcast1to16f(&ray16.tfar[rayIndex]);

	  const unsigned int leaf_mask = BVH4I_LEAF_MASK;

	  while (1)
	    {

	      NodeRef curNode = stack_node[sindex-1];
	      sindex--;

	      traverse_single_intersect<ENABLE_COMPRESSED_BVH4I_NODES>(curNode,
								      sindex,
								      rdir_xyz,
								      org_rdir_xyz,
								      min_dist_xyz,
								      max_dist_xyz,
								      stack_node,
								      stack_dist,
								      nodes,
								      leaf_mask);
		   


	      /* return if stack is empty */
	      if (unlikely(curNode == BVH4i::invalidNode)) break;

	      STAT3(normal.trav_leaves,1,1,1);
	      STAT3(normal.trav_prims,4,4,4);

	      /* intersect one ray against four triangles */

	      //////////////////////////////////////////////////////////////////////////////////////////////////
	      // PING;
	      const bool hit = LeafIntersector::intersect(curNode,
							  rayIndex,
							  dir_xyz,
							  org_xyz,
							  min_dist_xyz,
							  max_dist_xyz,
							  ray16,
							  accel,
							  (Scene*)bvh->geometry);
									   
	      if (hit)
		compactStack(stack_node,stack_dist,sindex,max_dist_xyz);

	      // ------------------------
	    }	  
	}
    }

    template<typename LeafIntersector,bool ENABLE_COMPRESSED_BVH4I_NODES>    
    void BVH4iIntersector16Test<LeafIntersector,ENABLE_COMPRESSED_BVH4I_NODES>::occluded(mic_i* valid_i, BVH4i* bvh, Ray16& ray16)
    {
      /* near and node stack */
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

      /* setup */
      const mic_m m_valid = *(mic_i*)valid_i != mic_i(0);
      const mic3f rdir16  = rcp_safe(ray16.dir);
      mic_m terminated    = !m_valid;
      const mic_f inf     = mic_f(pos_inf);
      const mic_f zero    = mic_f::zero();

      const Node      * __restrict__ nodes = (Node     *)bvh->nodePtr();
      const Triangle1 * __restrict__ accel = (Triangle1*)bvh->triPtr();

      stack_node[0] = BVH4i::invalidNode;

      long rayIndex = -1;
      while((rayIndex = bitscan64(rayIndex,toInt(m_valid))) != BITSCAN_NO_BIT_SET_64)	    
        {
	  stack_node[1] = bvh->root;
	  size_t sindex = 2;

	  const mic_f org_xyz      = loadAOS4to16f(rayIndex,ray16.org.x,ray16.org.y,ray16.org.z);
	  const mic_f dir_xyz      = loadAOS4to16f(rayIndex,ray16.dir.x,ray16.dir.y,ray16.dir.z);
	  const mic_f rdir_xyz     = loadAOS4to16f(rayIndex,rdir16.x,rdir16.y,rdir16.z);
	  const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
	  const mic_f min_dist_xyz = broadcast1to16f(&ray16.tnear[rayIndex]);
	  const mic_f max_dist_xyz = broadcast1to16f(&ray16.tfar[rayIndex]);
	  const mic_i v_invalidNode(BVH4i::invalidNode);
	  const unsigned int leaf_mask = BVH4I_LEAF_MASK;
	  const mic_i i_leaf_mask(BVH4I_LEAF_MASK);
	  while (1)
	    {
	      NodeRef curNode = stack_node[sindex-1];
	      sindex--;

	      const mic_m m7777 = 0x7777; 
	      const mic_m m_rdir0 = lt(m7777,rdir_xyz,mic_f::zero());
	      const mic_m m_rdir1 = ge(m7777,rdir_xyz,mic_f::zero());

	      while (1) 
		{
		  /* test if this is a leaf node */
		  if (unlikely(curNode.isLeaf(leaf_mask))) break;
		  STAT3(shadow.trav_nodes,1,1,1);

		  const BVH4i::Node* __restrict__ const node = curNode.node(nodes);

		  mic_f tLowerXYZ = select(m7777,rdir_xyz,min_dist_xyz); 
		  mic_f tUpperXYZ = select(m7777,rdir_xyz,max_dist_xyz);
		  mic_m hitm = ~m7777; 

		  const float* __restrict const plower = (float*)node->lower;
		  const float* __restrict const pupper = (float*)node->upper;

		  prefetch<PFHINT_L1>((char*)node + 0);
		  prefetch<PFHINT_L1>((char*)node + 64);
	    
		  /* intersect single ray with 4 bounding boxes */

		  tLowerXYZ = mask_msub(m_rdir1,tLowerXYZ,load16f(plower),org_rdir_xyz);
		  tUpperXYZ = mask_msub(m_rdir0,tUpperXYZ,load16f(plower),org_rdir_xyz);

		  tLowerXYZ = mask_msub(m_rdir0,tLowerXYZ,load16f(pupper),org_rdir_xyz);
		  tUpperXYZ = mask_msub(m_rdir1,tUpperXYZ,load16f(pupper),org_rdir_xyz);

		  const mic_f tLower = tLowerXYZ;
		  const mic_f tUpper = tUpperXYZ;


		  sindex--;
		  curNode = stack_node[sindex]; // early pop of next node

		  const mic_f tNear = vreduce_max4(tLower);
		  const mic_f tFar  = vreduce_min4(tUpper);  
		  hitm = le(hitm,tNear,tFar);
		  



		  STAT3(shadow.trav_hit_boxes[countbits(hitm)],1,1,1);


		  /* if no child is hit, continue with early popped child */
		  const mic_i plower_node = load16i((int*)node);

		  if (unlikely(none(hitm))) continue;
		  sindex++;
        
		  const unsigned long hiti = toInt(hitm);
		  const unsigned long pos_first = bitscan64(hiti);
		  const unsigned long num_hitm = countbits(hiti); 
        
		  /* if a single child is hit, continue with that child */
		  curNode = ((unsigned int *)node)[pos_first];
		  if (likely(num_hitm == 1)) continue;

		  /* if two children are hit, push in correct order */
		  const unsigned long pos_second = bitscan64(pos_first,hiti);
		  if (likely(num_hitm == 2))
		    {
		      const unsigned int dist_first  = ((unsigned int*)&tNear)[pos_first];
		      const unsigned int dist_second = ((unsigned int*)&tNear)[pos_second];
		      const unsigned int node_first  = curNode;
		      const unsigned int node_second = ((unsigned int*)node)[pos_second];
          
		      if (dist_first <= dist_second)
			{
			  stack_node[sindex] = node_second;
			  sindex++;
			  assert(sindex < 3*BVH4i::maxDepth+1);
			  continue;
			}
		      else
			{
			  stack_node[sindex] = curNode;
			  curNode = node_second;
			  sindex++;
			  assert(sindex < 3*BVH4i::maxDepth+1);
			  continue;
			}
		    }

		  /* continue with closest child and push all others */
		  // const mic_m m_leaf = test(hitm,plower_node,i_leaf_mask);
		  // hitm ^= m_leaf;
		  // const unsigned int num_nodes  = countbits(hitm);
		  // const unsigned int num_leaves = countbits(m_leaf);

		  // compactustore16i(hitm,&stack_node[sindex],plower_node);
		  // sindex += num_nodes;
		  // compactustore16i(m_leaf,&stack_node[sindex],plower_node);
		  // sindex += num_leaves;
		  // sindex -= 1;
		  // curNode = stack_node[sindex];

		  const mic_f tNear_pos = select(hitm,tNear,inf);
		  const mic_f min_dist = set_min_lanes(tNear_pos);
		  const unsigned int old_sindex = sindex;
		  sindex += countbits(hiti) - 1;

		  assert(sindex < 3*BVH4i::maxDepth+1);
        
		  const mic_m closest_child = eq(hitm,min_dist,tNear);
		  const unsigned long closest_child_pos = bitscan64(closest_child);
		  const mic_m m_pos = andn(hitm,andn(closest_child,(mic_m)((unsigned int)closest_child - 1)));
		  curNode = ((unsigned int*)node)[closest_child_pos];
		  compactustore16i(m_pos,&stack_node[old_sindex],plower_node);
		}

	      /* return if stack is empty */
	      if (unlikely(curNode == BVH4i::invalidNode)) break;

	      STAT3(shadow.trav_leaves,1,1,1);
	      STAT3(shadow.trav_prims,4,4,4);

	      /* intersect one ray against four triangles */

	      //////////////////////////////////////////////////////////////////////////////////////////////////

	      const bool hit = LeafIntersector::occluded(curNode,
							 rayIndex,
							 dir_xyz,
							 org_xyz,
							 min_dist_xyz,
							 max_dist_xyz,
							 ray16,
							 terminated,
							 accel,
							 (Scene*)bvh->geometry);

	      if (unlikely(hit)) break;
	      //////////////////////////////////////////////////////////////////////////////////////////////////

	    }


	  if (unlikely(all(toMask(terminated)))) break;
	}


      store16i(m_valid & toMask(terminated),&ray16.geomID,0);
    }

    typedef BVH4iIntersector16Test< Triangle1LeafIntersector  < true >, false  > Triangle1Intersector16TestMoellerFilter;
    typedef BVH4iIntersector16Test< Triangle1LeafIntersector  < false >, false > Triangle1Intersector16TestMoellerNoFilter;
    typedef BVH4iIntersector16Test< Triangle1mcLeafIntersector< true >, true  > Triangle1mcIntersector16TestMoellerFilter;
    typedef BVH4iIntersector16Test< Triangle1mcLeafIntersector< false >, true > Triangle1mcIntersector16TestMoellerNoFilter;


    DEFINE_INTERSECTOR16    (BVH4iTriangle1Intersector16TestMoeller          , Triangle1Intersector16TestMoellerFilter);
    DEFINE_INTERSECTOR16    (BVH4iTriangle1Intersector16TestMoellerNoFilter  , Triangle1Intersector16TestMoellerNoFilter);
    DEFINE_INTERSECTOR16    (BVH4iTriangle1mcIntersector16TestMoeller        , Triangle1mcIntersector16TestMoellerFilter);
    DEFINE_INTERSECTOR16    (BVH4iTriangle1mcIntersector16TestMoellerNoFilter, Triangle1mcIntersector16TestMoellerNoFilter);


  }
}


/*

      __aligned(64) mic_f   stack_dist[3*BVH4i::maxDepth+1];
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

      const mic_m m_valid = *(mic_i*)valid_i != mic_i(0);
      const unsigned int numValidRays = countbits(m_valid);
      __aligned(64) Vec3fa dir[16];
      __aligned(64) Vec3fa org[16];

      convertSOA4ftoAOS4f(ray16.org.x,ray16.org.y,ray16.org.z,ray16.tnear,org);
      convertSOA4ftoAOS4f(ray16.dir.x,ray16.dir.y,ray16.dir.z,ray16.tfar,dir);

      const mic_f inf        = mic_f(pos_inf);
      const mic_f zero       = mic_f::zero();

      store16f(stack_dist,inf);

      stack_node[0] = BVH4i::invalidNode;

      const Node      * __restrict__ nodes = (Node     *)bvh->nodePtr();
      const Triangle1 * __restrict__ accel = (Triangle1*)bvh->triPtr();
      for (size_t rayIndex4 = 0;rayIndex4<numValidRays;rayIndex4+=4)
        {
	  stack_node[1] = bvh->root;
	  size_t sindex = 2;
	  const mic_f dir4_xyz  = load16f(&dir[rayIndex4]);
	  const mic_f org4_xyz  = load16f(&org[rayIndex4]);
	  const mic_f rdir4_xyz = rcp_safe(dir4_xyz);
	  const mic_f min_dist4 = swDDDD(org4_xyz);
	  mic_f       max_dist4 = swDDDD(dir4_xyz);
	  const mic_f org4_rdir4_xyz = org4_xyz * rdir4_xyz;

	  const unsigned int leaf_mask = BVH4I_LEAF_MASK;

	  while (1)
	    {

	      NodeRef curNode = stack_node[sindex-1];
	      mic_f curDist   = stack_dist[sindex-1];
	      sindex--;
	      const mic_m m_stackDist = max_dist4 > curDist;
	      if (unlikely(curNode == BVH4i::invalidNode))  break;
        
	      if (unlikely(none(m_stackDist))) {continue;}

	      const mic_m m7777 = 0x7777; 

	      const mic_f org_tLowerXYZ = select(m7777,rdir4_xyz,min_dist4); 
	      const mic_f org_tUpperXYZ = select(m7777,rdir4_xyz,max_dist4);
    
	      while (1)
		{
		  if (unlikely(curNode.isLeaf(leaf_mask))) break;
          
		  //STAT3(normal.trav_nodes,1,popcnt(ray_tfar > curDist),16);
		  const BVH4i::Node* __restrict__ const node = curNode.node(nodes);


		  prefetch<PFHINT_L1>((mic_f*)node + 0);           
		  prefetch<PFHINT_L1>((mic_f*)node + 1); 

		  sindex--;

		  curNode = stack_node[sindex]; 	  
		  curDist = stack_dist[sindex];

#pragma unroll(4)
		  for (unsigned int i=0; i<4; i++)
		    {
		      BVH4i::NodeRef child = node->lower[i].child;
		      
		      const mic_f lower = broadcast4to16f(&node->lower[i]);
		      const mic_f upper = broadcast4to16f(&node->upper[i]);

		      mic_f tLowerXYZ = org_tLowerXYZ;		      
		      mic_f tUpperXYZ = org_tUpperXYZ;

		      tLowerXYZ = mask_msub(m7777,tLowerXYZ,lower,org4_rdir4_xyz);
		      tUpperXYZ = mask_msub(m7777,tUpperXYZ,upper,org4_rdir4_xyz);

		      if (unlikely(i >=2 && child == BVH4i::invalidNode)) break;
	    
		      const mic_f tLower = min(tLowerXYZ,tUpperXYZ);
		      const mic_f tUpper = max(tLowerXYZ,tUpperXYZ);

		      const mic_f tNear = vreduce_max4(tLower);
		      const mic_f tFar  = vreduce_min4(tUpper);  

		      const mic_m hitm = le(0x8888,tNear,tFar);


		      const mic_f childDist = select(hitm,tNear,inf);
		      const mic_m m_child_dist = lt(childDist,curDist);

		      if (likely(any(hitm)))
			{
			  sindex++;
			  if (any(m_child_dist))
			    {
			      stack_node[sindex-1] = curNode;
			      stack_dist[sindex-1] = curDist; 
			      curDist = childDist;
			      curNode = child;
			    }              
			  else 
			    {
			      stack_node[sindex-1] = child;
			      stack_dist[sindex-1] = childDist; 
			    }
			}	      

		    }
		}		   


	      if (unlikely(curNode == BVH4i::invalidNode)) break;

	      STAT3(normal.trav_leaves,1,1,1);
	      STAT3(normal.trav_prims,4,4,4);


	      //////////////////////////////////////////////////////////////////////////////////////////////////

	      const mic_m m_valid_leaf = 0xffff;
	      STAT3(normal.trav_leaves,1,popcnt(m_valid_leaf),16);
 
	      LeafIntersector::intersect16(curNode,m_valid_leaf,ray16.dir,ray16.org,ray16,accel,(Scene*)bvh->geometry);

	      //ray_tfar = select(m_valid_leaf,ray.tfar,ray_tfar);


	      // ------------------------
	    }	  
	}

*/
