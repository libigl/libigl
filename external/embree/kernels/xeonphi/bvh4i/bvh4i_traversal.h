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

#pragma once

#include "bvh4i.h"

namespace embree
{
  // ====================================================================================================================
  // ====================================================================================================================
  // ====================================================================================================================

  template<bool DECOMPRESS_NODE>
  __forceinline void traverse_single_intersect(BVH4i::NodeRef &curNode,
					       size_t &sindex,
					       const mic_f &rdir_xyz,
					       const mic_f &org_rdir_xyz,
					       const mic_f &min_dist_xyz,
					       const mic_f &max_dist_xyz,
					       BVH4i::NodeRef *__restrict__ const stack_node,
					       float   *__restrict__ const stack_dist,
					       const BVH4i::Node      * __restrict__ const nodes,
					       const unsigned int leaf_mask)
  {
    const mic_m m7777 = 0x7777; 
    const mic_m m_rdir0 = lt(m7777,rdir_xyz,mic_f::zero());
    const mic_m m_rdir1 = ge(m7777,rdir_xyz,mic_f::zero());

    while (1) 
      {
	/* test if this is a leaf node */
	if (unlikely(curNode.isLeaf(leaf_mask))) break;
	STAT3(normal.trav_nodes,1,1,1);

	const BVH4i::Node* __restrict__ const node = curNode.node(nodes);

	mic_f tLowerXYZ = select(m7777,rdir_xyz,min_dist_xyz); 
	mic_f tUpperXYZ = select(m7777,rdir_xyz,max_dist_xyz);
	mic_m hitm = ~m7777; 

	if (!DECOMPRESS_NODE)
	  {
	    const float* __restrict const plower = (float*)node->lower;
	    const float* __restrict const pupper = (float*)node->upper;

	    prefetch<PFHINT_L1>((char*)node + 0);
	    prefetch<PFHINT_L1>((char*)node + 64);
	    
	    /* intersect single ray with 4 bounding boxes */

	    tLowerXYZ = mask_msub(m_rdir1,tLowerXYZ,load16f(plower),org_rdir_xyz);
	    tUpperXYZ = mask_msub(m_rdir0,tUpperXYZ,load16f(plower),org_rdir_xyz);

	    tLowerXYZ = mask_msub(m_rdir0,tLowerXYZ,load16f(pupper),org_rdir_xyz);
	    tUpperXYZ = mask_msub(m_rdir1,tUpperXYZ,load16f(pupper),org_rdir_xyz);
	  }
	else
	  {
	    BVH4i::QuantizedNode* __restrict__ const compressed_node = (BVH4i::QuantizedNode*)node;
	    prefetch<PFHINT_L1>((char*)node + 0);
		  
	    //DBG_PRINT(*compressed_node);
	    const mic_f startXYZ = compressed_node->decompress_startXYZ();
	    const mic_f diffXYZ  = compressed_node->decompress_diffXYZ();
	    const mic_f clower   = compressed_node->decompress_lowerXYZ(startXYZ,diffXYZ);
	    const mic_f cupper   = compressed_node->decompress_upperXYZ(startXYZ,diffXYZ);

	    tLowerXYZ = mask_msub(m_rdir1,tLowerXYZ,clower,org_rdir_xyz);
	    tUpperXYZ = mask_msub(m_rdir0,tUpperXYZ,clower,org_rdir_xyz);

	    tLowerXYZ = mask_msub(m_rdir0,tLowerXYZ,cupper,org_rdir_xyz);
	    tUpperXYZ = mask_msub(m_rdir1,tUpperXYZ,cupper,org_rdir_xyz);	    
	  }

	const mic_f tLower = tLowerXYZ;
	const mic_f tUpper = tUpperXYZ;

	/* early pop of next node */
	sindex--;
	curNode = stack_node[sindex];

#ifdef RTCORE_STAT_COUNTERS
	if (!curNode.isLeaf(leaf_mask))
	  STAT3(normal.trav_stack_nodes,1,1,1);
#endif


	const mic_f tNear = vreduce_max4(tLower);
	const mic_f tFar  = vreduce_min4(tUpper);  
	hitm = le(hitm,tNear,tFar);
		  
	const mic_f tNear_pos = select(hitm,tNear,inf);

	STAT3(normal.trav_hit_boxes[countbits(hitm)],1,1,1);

	const mic_i plower_node = load16i((int*)node);


	/* if no child is hit, continue with early popped child */
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
		((unsigned int*)stack_dist)[sindex] = dist_second;                      
		sindex++;
		assert(sindex < 3*BVH4i::maxDepth+1);
		continue;
	      }
	    else
	      {
		stack_node[sindex] = curNode;
		((unsigned int*)stack_dist)[sindex] = dist_first;
		curNode = node_second;
		sindex++;
		assert(sindex < 3*BVH4i::maxDepth+1);
		continue;
	      }
	  }

	/* continue with closest child and push all others */


	const mic_f min_dist = set_min_lanes(tNear_pos);
	const unsigned int old_sindex = sindex;
	sindex += countbits(hiti) - 1;
	assert(sindex < 3*BVH4i::maxDepth+1);
        
	const BVH4i::Node* __restrict__ const next = curNode.node(nodes);
	prefetch<PFHINT_L1>((char*)next + 0*64);
	prefetch<PFHINT_L1>((char*)next + 1*64);

	const mic_m closest_child = eq(hitm,min_dist,tNear);
	const unsigned long closest_child_pos = bitscan64(closest_child);
	const mic_m m_pos = andn(hitm,andn(closest_child,(mic_m)((unsigned int)closest_child - 1)));
	curNode = ((unsigned int*)node)[closest_child_pos];
	compactustore16f(m_pos,&stack_dist[old_sindex],tNear);
	compactustore16i(m_pos,&stack_node[old_sindex],plower_node);

	if (unlikely(((unsigned int*)stack_dist)[sindex-2] < ((unsigned int*)stack_dist)[sindex-1]))
	  {
	    std::swap(((unsigned int*)stack_dist)[sindex-1],((unsigned int*)stack_dist)[sindex-2]);
	    std::swap(((unsigned int*)stack_node)[sindex-1],((unsigned int*)stack_node)[sindex-2]);
	  }
      }

  }

  template<bool DECOMPRESS_NODE>
  __forceinline void traverse_single_occluded(BVH4i::NodeRef &curNode,
					      size_t &sindex,
					      const mic_f &rdir_xyz,
					      const mic_f &org_rdir_xyz,
					      const mic_f &min_dist_xyz,
					      const mic_f &max_dist_xyz,
					      BVH4i::NodeRef *__restrict__ const stack_node,
					      const BVH4i::Node      * __restrict__ const nodes,
					      const unsigned int leaf_mask)
  {
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

	if (!DECOMPRESS_NODE)
	  {
	    const float* __restrict const plower = (float*)node->lower;
	    const float* __restrict const pupper = (float*)node->upper;

	    prefetch<PFHINT_L1>((char*)node + 0);
	    prefetch<PFHINT_L1>((char*)node + 64);
	    
	    /* intersect single ray with 4 bounding boxes */

	    tLowerXYZ = mask_msub(m_rdir1,tLowerXYZ,load16f(plower),org_rdir_xyz);
	    tUpperXYZ = mask_msub(m_rdir0,tUpperXYZ,load16f(plower),org_rdir_xyz);

	    tLowerXYZ = mask_msub(m_rdir0,tLowerXYZ,load16f(pupper),org_rdir_xyz);
	    tUpperXYZ = mask_msub(m_rdir1,tUpperXYZ,load16f(pupper),org_rdir_xyz);
	  }
	else
	  {
	    BVH4i::QuantizedNode* __restrict__ const compressed_node = (BVH4i::QuantizedNode*)node;
	    prefetch<PFHINT_L1>((char*)node + 0);
		  
	    const mic_f startXYZ = compressed_node->decompress_startXYZ();
	    const mic_f diffXYZ  = compressed_node->decompress_diffXYZ();
	    const mic_f clower   = compressed_node->decompress_lowerXYZ(startXYZ,diffXYZ);
	    const mic_f cupper   = compressed_node->decompress_upperXYZ(startXYZ,diffXYZ);

	    tLowerXYZ = mask_msub(m_rdir1,tLowerXYZ,clower,org_rdir_xyz);
	    tUpperXYZ = mask_msub(m_rdir0,tUpperXYZ,clower,org_rdir_xyz);

	    tLowerXYZ = mask_msub(m_rdir0,tLowerXYZ,cupper,org_rdir_xyz);
	    tUpperXYZ = mask_msub(m_rdir1,tUpperXYZ,cupper,org_rdir_xyz);	    
	  }

	const mic_f tLower = tLowerXYZ;
	const mic_f tUpper = tUpperXYZ;


	sindex--;
	curNode = stack_node[sindex]; // early pop of next node

#ifdef RTCORE_STAT_COUNTERS
	if (!curNode.isLeaf(leaf_mask))
	  STAT3(shadow.trav_stack_nodes,1,1,1);
#endif

	const mic_f tNear = vreduce_max4(tLower);
	const mic_f tFar  = vreduce_min4(tUpper);  
	hitm = le(hitm,tNear,tFar);
		  

	const mic_f tNear_pos = select(hitm,tNear,inf);


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

  }

  __forceinline void compactStack(BVH4i::NodeRef *__restrict__ const stack_node,
				  float   *__restrict__ const stack_dist,
				  size_t &sindex,
				  const mic_f &max_dist_xyz)
  {
    if (likely(sindex >= 2))
      {
	if (likely(sindex < 16))
	  {
	    const unsigned int m_num_stack = mic_m::shift1[sindex] - 1;
	    const mic_m m_num_stack_low  = toMask(m_num_stack);
	    const mic_f snear_low  = load16f(stack_dist + 0);
	    const mic_i snode_low  = load16i((int*)stack_node + 0);
	    const mic_m m_stack_compact_low  = le(m_num_stack_low,snear_low,max_dist_xyz) | (mic_m)1;
	    compactustore16f_low(m_stack_compact_low,stack_dist + 0,snear_low);
	    compactustore16i_low(m_stack_compact_low,(int*)stack_node + 0,snode_low);
	    sindex = countbits(m_stack_compact_low);
	    assert(sindex < 16);
	  }
	else if (likely(sindex < 32))
	  {
	    const mic_m m_num_stack_high = toMask(mic_m::shift1[sindex-16] - 1); 
	    const mic_f snear_low  = load16f(stack_dist + 0);
	    const mic_f snear_high = load16f(stack_dist + 16);
	    const mic_i snode_low  = load16i((int*)stack_node + 0);
	    const mic_i snode_high = load16i((int*)stack_node + 16);
	    const mic_m m_stack_compact_low  = le(snear_low,max_dist_xyz) | (mic_m)1;
	    const mic_m m_stack_compact_high = le(m_num_stack_high,snear_high,max_dist_xyz);
	    compactustore16f(m_stack_compact_low,      stack_dist + 0,snear_low);
	    compactustore16i(m_stack_compact_low,(int*)stack_node + 0,snode_low);
	    compactustore16f(m_stack_compact_high,      stack_dist + countbits(m_stack_compact_low),snear_high);
	    compactustore16i(m_stack_compact_high,(int*)stack_node + countbits(m_stack_compact_low),snode_high);
	    assert ((unsigned int )m_num_stack_high == ((mic_m::shift1[sindex] - 1) >> 16));

	    sindex = countbits(m_stack_compact_low) + countbits(m_stack_compact_high);
	    assert(sindex < 32);
	  }
	else
	  {
	    const mic_m m_num_stack_32 = toMask(mic_m::shift1[sindex-32] - 1); 

	    const mic_f snear_0  = load16f(stack_dist + 0);
	    const mic_f snear_16 = load16f(stack_dist + 16);
	    const mic_f snear_32 = load16f(stack_dist + 32);
	    const mic_i snode_0  = load16i((int*)stack_node + 0);
	    const mic_i snode_16 = load16i((int*)stack_node + 16);
	    const mic_i snode_32 = load16i((int*)stack_node + 32);
	    const mic_m m_stack_compact_0  = le(               snear_0 ,max_dist_xyz) | (mic_m)1;
	    const mic_m m_stack_compact_16 = le(               snear_16,max_dist_xyz);
	    const mic_m m_stack_compact_32 = le(m_num_stack_32,snear_32,max_dist_xyz);

	    sindex = 0;
	    compactustore16f(m_stack_compact_0,      stack_dist + sindex,snear_0);
	    compactustore16i(m_stack_compact_0,(int*)stack_node + sindex,snode_0);
	    sindex += countbits(m_stack_compact_0);
	    compactustore16f(m_stack_compact_16,      stack_dist + sindex,snear_16);
	    compactustore16i(m_stack_compact_16,(int*)stack_node + sindex,snode_16);
	    sindex += countbits(m_stack_compact_16);
	    compactustore16f(m_stack_compact_32,      stack_dist + sindex,snear_32);
	    compactustore16i(m_stack_compact_32,(int*)stack_node + sindex,snode_32);
	    sindex += countbits(m_stack_compact_32);

	    assert(sindex < 48);		  
	  }
      }
  }

  template<bool DECOMPRESS_NODE>
  __forceinline void traverse_chunk_intersect(BVH4i::NodeRef &curNode,
					      mic_f &curDist,
					      const mic3f &rdir,
					      const mic3f &org_rdir,
					      const mic_f &ray_tnear,
					      const mic_f &ray_tfar,
					      BVH4i::NodeRef *__restrict__ &sptr_node,
					      mic_f *__restrict__ &sptr_dist,
					      const BVH4i::Node      * __restrict__ const nodes,
					      const unsigned int leaf_mask)
  {

    while (1)
      {
	/* test if this is a leaf node */
	if (unlikely(curNode.isLeaf(leaf_mask))) break;
          
	STAT3(normal.trav_nodes,1,popcnt(ray_tfar > curDist),16);
	const BVH4i::Node* __restrict__ const node = curNode.node(nodes);


	/* pop of next node */
	sptr_node--;
	sptr_dist--;
	curNode = *sptr_node; 	  
	curDist = *sptr_dist;

	prefetch<PFHINT_L1>((mic_f*)node + 0);           
	prefetch<PFHINT_L1>((mic_f*)node + 1); 

#pragma unroll(4)
	for (unsigned int i=0; i<4; i++)
          {
	    BVH4i::NodeRef child;
	    mic_f lclipMinX,lclipMinY,lclipMinZ;
	    mic_f lclipMaxX,lclipMaxY,lclipMaxZ;

	    if (!DECOMPRESS_NODE)
	      {
		child = node->lower[i].child;

		lclipMinX = msub(node->lower[i].x,rdir.x,org_rdir.x);
		lclipMinY = msub(node->lower[i].y,rdir.y,org_rdir.y);
		lclipMinZ = msub(node->lower[i].z,rdir.z,org_rdir.z);
		lclipMaxX = msub(node->upper[i].x,rdir.x,org_rdir.x);
		lclipMaxY = msub(node->upper[i].y,rdir.y,org_rdir.y);
		lclipMaxZ = msub(node->upper[i].z,rdir.z,org_rdir.z);
	      }
	    else
	      {
		BVH4i::QuantizedNode* __restrict__ const compressed_node = (BVH4i::QuantizedNode*)node;
		child = compressed_node->child(i);

		const mic_f startXYZ = compressed_node->decompress_startXYZ();
		const mic_f diffXYZ  = compressed_node->decompress_diffXYZ();
		const mic_f clower   = compressed_node->decompress_lowerXYZ(startXYZ,diffXYZ);
		const mic_f cupper   = compressed_node->decompress_upperXYZ(startXYZ,diffXYZ);

		lclipMinX = msub(mic_f(clower[4*i+0]),rdir.x,org_rdir.x);
		lclipMinY = msub(mic_f(clower[4*i+1]),rdir.y,org_rdir.y);
		lclipMinZ = msub(mic_f(clower[4*i+2]),rdir.z,org_rdir.z);
		lclipMaxX = msub(mic_f(cupper[4*i+0]),rdir.x,org_rdir.x);
		lclipMaxY = msub(mic_f(cupper[4*i+1]),rdir.y,org_rdir.y);
		lclipMaxZ = msub(mic_f(cupper[4*i+2]),rdir.z,org_rdir.z);		
	      }

	    if (unlikely(i >=2 && child == BVH4i::invalidNode)) break;
	    
            const mic_f lnearP = max(max(min(lclipMinX, lclipMaxX), min(lclipMinY, lclipMaxY)), min(lclipMinZ, lclipMaxZ));
            const mic_f lfarP  = min(min(max(lclipMinX, lclipMaxX), max(lclipMinY, lclipMaxY)), max(lclipMinZ, lclipMaxZ));
            const mic_m lhit   = le(max(lnearP,ray_tnear),min(lfarP,ray_tfar));   
	    const mic_f childDist = select(lhit,lnearP,inf);
            const mic_m m_child_dist = lt(childDist,curDist);
            /* if we hit the child we choose to continue with that child if it 
               is closer than the current next child, or we push it onto the stack */


            if (likely(any(lhit)))
	      {
		sptr_node++;
		sptr_dist++;
              
		/* push cur node onto stack and continue with hit child */
		if (any(m_child_dist))
		  {
		    *(sptr_node-1) = curNode;
		    *(sptr_dist-1) = curDist; 
		    curDist = childDist;
		    curNode = child;
		  }              
		/* push hit child onto stack*/
		else 
		  {
		    *(sptr_node-1) = child;
		    *(sptr_dist-1) = childDist; 
		  }
	      }	      

          }
      }

  }


  template<bool DECOMPRESS_NODE>
  __forceinline void traverse_chunk_occluded(BVH4i::NodeRef &curNode,
					     mic_f &curDist,
					     const mic3f &rdir,
					     const mic3f &org_rdir,
					     const mic_f &ray_tnear,
					     const mic_f &ray_tfar,
					     const mic_m &m_active,
					     BVH4i::NodeRef *__restrict__ &sptr_node,
					     mic_f *__restrict__ &sptr_dist,
					     const BVH4i::Node      * __restrict__ const nodes,
					     const unsigned int leaf_mask)
  {
    while (1)
      {
	/* test if this is a leaf node */
	if (unlikely(curNode.isLeaf(leaf_mask))) break;
          
	STAT3(shadow.trav_nodes,1,popcnt(ray_tfar > curDist),16);
	const BVH4i::Node* __restrict__ const node = curNode.node(nodes);
          
	prefetch<PFHINT_L1>((mic_f*)node + 0); 
	prefetch<PFHINT_L1>((mic_f*)node + 1); 

	/* pop of next node */
	sptr_node--;
	sptr_dist--;
	curNode = *sptr_node; 
	curDist = *sptr_dist;
          	 
#pragma unroll(4)
	for (unsigned int i=0; i<4; i++)
          {
	    BVH4i::NodeRef child;
	    mic_f lclipMinX,lclipMinY,lclipMinZ;
	    mic_f lclipMaxX,lclipMaxY,lclipMaxZ;

	    if (!DECOMPRESS_NODE)
	      {
		child = node->lower[i].child;

		lclipMinX = msub(node->lower[i].x,rdir.x,org_rdir.x);
		lclipMinY = msub(node->lower[i].y,rdir.y,org_rdir.y);
		lclipMinZ = msub(node->lower[i].z,rdir.z,org_rdir.z);
		lclipMaxX = msub(node->upper[i].x,rdir.x,org_rdir.x);
		lclipMaxY = msub(node->upper[i].y,rdir.y,org_rdir.y);
		lclipMaxZ = msub(node->upper[i].z,rdir.z,org_rdir.z);
	      }
	    else
	      {
		BVH4i::QuantizedNode* __restrict__ const compressed_node = (BVH4i::QuantizedNode*)node;
		child = compressed_node->child(i);

		const mic_f startXYZ = compressed_node->decompress_startXYZ();
		const mic_f diffXYZ  = compressed_node->decompress_diffXYZ();
		const mic_f clower   = compressed_node->decompress_lowerXYZ(startXYZ,diffXYZ);
		const mic_f cupper   = compressed_node->decompress_upperXYZ(startXYZ,diffXYZ);

		lclipMinX = msub(mic_f(clower[4*i+0]),rdir.x,org_rdir.x);
		lclipMinY = msub(mic_f(clower[4*i+1]),rdir.y,org_rdir.y);
		lclipMinZ = msub(mic_f(clower[4*i+2]),rdir.z,org_rdir.z);
		lclipMaxX = msub(mic_f(cupper[4*i+0]),rdir.x,org_rdir.x);
		lclipMaxY = msub(mic_f(cupper[4*i+1]),rdir.y,org_rdir.y);
		lclipMaxZ = msub(mic_f(cupper[4*i+2]),rdir.z,org_rdir.z);		
	      }

	    if (unlikely(i >=2 && child == BVH4i::invalidNode)) break;

            const mic_f lnearP = max(max(min(lclipMinX, lclipMaxX), min(lclipMinY, lclipMaxY)), min(lclipMinZ, lclipMaxZ));
            const mic_f lfarP  = min(min(max(lclipMinX, lclipMaxX), max(lclipMinY, lclipMaxY)), max(lclipMinZ, lclipMaxZ));
            const mic_m lhit   = le(m_active,max(lnearP,ray_tnear),min(lfarP,ray_tfar));      
	    const mic_f childDist = select(lhit,lnearP,inf);
            const mic_m m_child_dist = childDist < curDist;
            
            /* if we hit the child we choose to continue with that child if it 
               is closer than the current next child, or we push it onto the stack */
            if (likely(any(lhit)))
	      {
		sptr_node++;
		sptr_dist++;
              
		/* push cur node onto stack and continue with hit child */
		if (any(m_child_dist))
		  {
		    *(sptr_node-1) = curNode; 
		    *(sptr_dist-1) = curDist; 
		    curDist = childDist;
		    curNode = child;
		  }
              
		/* push hit child onto stack*/
		else {
		  *(sptr_node-1) = child;
		  *(sptr_dist-1) = childDist; 
		}
	      }	      
          }
      }
    
  }


  /* BVH4i::QuantizedNode* __restrict__ const node = (BVH4i::QuantizedNode*)curNode.node(nodes); */
  /* prefetch<PFHINT_L1>((char*)node + 0); */

  /* const float* __restrict const plower = (float*)node; */
		  
  /* const mic_f startXYZ = node->decompress_startXYZ(); */
  /* const mic_f diffXYZ  = node->decompress_diffXYZ(); */
  /* const mic_f lower   = node->decompress_lowerXYZ(startXYZ,diffXYZ); */
  /* const mic_f upper   = node->decompress_upperXYZ(startXYZ,diffXYZ); */

  /* mic_f tLowerXYZ = lower * rdir_xyz - org_rdir_xyz; */
  /* mic_f tUpperXYZ = upper * rdir_xyz - org_rdir_xyz; */


  /* BVH4i::QuantizedNode* __restrict__ const node = (BVH4i::QuantizedNode*)curNode.node(nodes); */
  /* prefetch<PFHINT_L1>((char*)node + 0); */

  /* const float* __restrict const plower = (float*)node; */
		  
  /* const mic_f startXYZ = node->decompress_startXYZ(); */
  /* const mic_f diffXYZ  = node->decompress_diffXYZ(); */
  /* const mic_f lower   = node->decompress_lowerXYZ(startXYZ,diffXYZ); */
  /* const mic_f upper   = node->decompress_upperXYZ(startXYZ,diffXYZ); */

  /* const mic_f tLowerXYZ = lower * rdir_xyz - org_rdir_xyz; */
  /* const mic_f tUpperXYZ = upper * rdir_xyz - org_rdir_xyz; */
  
};
