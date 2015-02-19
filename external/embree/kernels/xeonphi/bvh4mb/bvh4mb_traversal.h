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

#include "bvh4mb.h"
#include "bvh4i/bvh4i_traversal.h"

namespace embree
{
  __forceinline void traverse_single_intersect(BVH4i::NodeRef &curNode,
					       size_t &sindex,
					       const mic_f &rdir_xyz,
					       const mic_f &org_rdir_xyz,
					       const mic_f &min_dist_xyz,
					       const mic_f &max_dist_xyz,
					       const mic_f &time,
					       BVH4i::NodeRef *__restrict__ const stack_node,
					       float   *__restrict__ const stack_dist,
					       const BVH4i::Node      * __restrict__ const nodes,
					       const unsigned int leaf_mask)
  {
    const mic_f one_time = (mic_f::one() - time);
    const mic_m m7777 = 0x7777; 
    const mic_m m_rdir0 = lt(m7777,rdir_xyz,mic_f::zero());
    const mic_m m_rdir1 = ge(m7777,rdir_xyz,mic_f::zero());

    while (1) 
      {
	/* test if this is a leaf node */
	if (unlikely(curNode.isLeaf(leaf_mask))) break;
        
	const BVH4i::Node* __restrict__ const node = curNode.node(nodes);
	const float* __restrict const plower = (float*)node->lower;
	const float* __restrict const pupper = (float*)node->upper;

	prefetch<PFHINT_L1>((char*)node + 0*64);
	prefetch<PFHINT_L1>((char*)node + 1*64);
	prefetch<PFHINT_L1>((char*)node + 2*64);
	prefetch<PFHINT_L1>((char*)node + 3*64);

	const BVH4mb::Node* __restrict__ const nodeMB = (BVH4mb::Node*)node;
	const mic_f lower = one_time  * load16f((float*)nodeMB->lower) + time * load16f((float*)nodeMB->lower_t1);
	const mic_f upper = one_time  * load16f((float*)nodeMB->upper) + time * load16f((float*)nodeMB->upper_t1);
		  
        
	/* intersect single ray with 4 bounding boxes */
	mic_f tLowerXYZ = select(m7777,rdir_xyz,min_dist_xyz);
	mic_f tUpperXYZ = select(m7777,rdir_xyz,max_dist_xyz);

	tLowerXYZ = mask_msub(m_rdir1,tLowerXYZ,lower,org_rdir_xyz);
	tUpperXYZ = mask_msub(m_rdir0,tUpperXYZ,lower,org_rdir_xyz);

	tLowerXYZ = mask_msub(m_rdir0,tLowerXYZ,upper,org_rdir_xyz);
	tUpperXYZ = mask_msub(m_rdir1,tUpperXYZ,upper,org_rdir_xyz);

	mic_m hitm = ~m7777; 
	const mic_f tLower = tLowerXYZ;
	const mic_f tUpper = tUpperXYZ;

	sindex--;

	curNode = stack_node[sindex]; // early pop of next node

	const BVH4i::Node* __restrict__ const next = curNode.node(nodes);
	prefetch<PFHINT_L2>((char*)next + 0);
	prefetch<PFHINT_L2>((char*)next + 64);

	const mic_f tNear = vreduce_max4(tLower);
	const mic_f tFar  = vreduce_min4(tUpper);  
	hitm = le(hitm,tNear,tFar);
	const mic_f tNear_pos = select(hitm,tNear,inf);


	/* if no child is hit, continue with early popped child */
	if (unlikely(none(hitm))) continue;
	sindex++;
        
	const unsigned long hiti = toInt(hitm);
	const unsigned long pos_first = bitscan64(hiti);
	const unsigned long num_hitm = countbits(hiti); 
        
	/* if a single child is hit, continue with that child */
	curNode = ((unsigned int *)plower)[pos_first];
	if (likely(num_hitm == 1)) continue;
        
	/* if two children are hit, push in correct order */
	const unsigned long pos_second = bitscan64(pos_first,hiti);
	if (likely(num_hitm == 2))
	  {
	    const unsigned int dist_first  = ((unsigned int*)&tNear)[pos_first];
	    const unsigned int dist_second = ((unsigned int*)&tNear)[pos_second];
	    const unsigned int node_first  = curNode;
	    const unsigned int node_second = ((unsigned int*)plower)[pos_second];
          
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
	const mic_i plower_node = load16i((int*)plower);
        
	const mic_m closest_child = eq(hitm,min_dist,tNear);
	const unsigned long closest_child_pos = bitscan64(closest_child);
	const mic_m m_pos = andn(hitm,andn(closest_child,(mic_m)((unsigned int)closest_child - 1)));
	curNode = ((unsigned int*)plower)[closest_child_pos];
	compactustore16f(m_pos,&stack_dist[old_sindex],tNear);
	compactustore16i(m_pos,&stack_node[old_sindex],plower_node);
      }
  }

  __forceinline void traverse_single_occluded(BVH4i::NodeRef &curNode,
					      size_t &sindex,
					      const mic_f &rdir_xyz,
					      const mic_f &org_rdir_xyz,
					      const mic_f &min_dist_xyz,
					      const mic_f &max_dist_xyz,
					      const mic_f &time,
					      BVH4i::NodeRef *__restrict__ const stack_node,
					      const BVH4i::Node      * __restrict__ const nodes,
					      const unsigned int leaf_mask)
  {
    const mic_f one_time = (mic_f::one() - time);
    const mic_m m7777 = 0x7777; 
    const mic_m m_rdir0 = lt(m7777,rdir_xyz,mic_f::zero());
    const mic_m m_rdir1 = ge(m7777,rdir_xyz,mic_f::zero());

    while (1) 
      {
	/* test if this is a leaf node */
	if (unlikely(curNode.isLeaf(leaf_mask))) break;
        
	const BVH4i::Node* __restrict__ const node = curNode.node(nodes);
	const float* __restrict const plower = (float*)node->lower;
	const float* __restrict const pupper = (float*)node->upper;

	prefetch<PFHINT_L1>((char*)node + 0*64);
	prefetch<PFHINT_L1>((char*)node + 1*64);
	prefetch<PFHINT_L1>((char*)node + 2*64);
	prefetch<PFHINT_L1>((char*)node + 3*64);

	const BVH4mb::Node* __restrict__ const nodeMB = (BVH4mb::Node*)node;
	const mic_f lower = one_time  * load16f((float*)nodeMB->lower) + time * load16f((float*)nodeMB->lower_t1);
	const mic_f upper = one_time  * load16f((float*)nodeMB->upper) + time * load16f((float*)nodeMB->upper_t1);
        
	/* intersect single ray with 4 bounding boxes */
	mic_f tLowerXYZ = select(m7777,rdir_xyz,min_dist_xyz);
	mic_f tUpperXYZ = select(m7777,rdir_xyz,max_dist_xyz);

	tLowerXYZ = mask_msub(m_rdir1,tLowerXYZ,lower,org_rdir_xyz);
	tUpperXYZ = mask_msub(m_rdir0,tUpperXYZ,lower,org_rdir_xyz);

	tLowerXYZ = mask_msub(m_rdir0,tLowerXYZ,upper,org_rdir_xyz);
	tUpperXYZ = mask_msub(m_rdir1,tUpperXYZ,upper,org_rdir_xyz);

	mic_m hitm = ~m7777; 
	const mic_f tLower = tLowerXYZ;
	const mic_f tUpper = tUpperXYZ;

	const BVH4i::Node* __restrict__ const next = curNode.node(nodes);
	prefetch<PFHINT_L2>((char*)next + 0);
	prefetch<PFHINT_L2>((char*)next + 64);

	sindex--;
	const mic_f tNear = vreduce_max4(tLower);
	const mic_f tFar  = vreduce_min4(tUpper);  
	hitm = le(hitm,tNear,tFar);
	const mic_f tNear_pos = select(hitm,tNear,inf);

	curNode = stack_node[sindex]; // early pop of next node

	/* if no child is hit, continue with early popped child */
	if (unlikely(none(hitm))) continue;
	sindex++;
        
	const unsigned long hiti = toInt(hitm);
	const unsigned long pos_first = bitscan64(hiti);
	const unsigned long num_hitm = countbits(hiti); 
        
	/* if a single child is hit, continue with that child */
	curNode = ((unsigned int *)plower)[pos_first];
	if (likely(num_hitm == 1)) continue;
        
	/* if two children are hit, push in correct order */
	const unsigned long pos_second = bitscan64(pos_first,hiti);
	if (likely(num_hitm == 2))
	  {
	    const unsigned int dist_first  = ((unsigned int*)&tNear)[pos_first];
	    const unsigned int dist_second = ((unsigned int*)&tNear)[pos_second];
	    const unsigned int node_first  = curNode;
	    const unsigned int node_second = ((unsigned int*)plower)[pos_second];
          
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
	const mic_i plower_node = load16i((int*)plower);
	curNode = ((unsigned int*)plower)[closest_child_pos];
	compactustore16i(m_pos,&stack_node[old_sindex],plower_node);
      }    
  }

  __forceinline void traverse_chunk_intersect(BVH4i::NodeRef &curNode,
					      mic_f &curDist,
					      const mic3f &rdir,
					      const mic3f &org_rdir,
					      const mic_f &ray_tnear,
					      const mic_f &ray_tfar,
					      const mic_f &time,
					      BVH4i::NodeRef *__restrict__ &sptr_node,
					      mic_f *__restrict__ &sptr_dist,
					      const BVH4i::Node      * __restrict__ const nodes,
					      const unsigned int leaf_mask)
  {


    const mic_f one_time = (mic_f::one() - time);
    while (1)
      {
	/* test if this is a leaf node */
	if (unlikely(curNode.isLeaf(leaf_mask))) break;
          
	STAT3(normal.trav_nodes,1,popcnt(ray_tfar > curDist),16);
	const BVH4i::Node* __restrict__ const node = curNode.node(nodes);

	const BVH4mb::Node* __restrict__ const nodeMB = (BVH4mb::Node*)node;

	/* pop of next node */
	sptr_node--;
	sptr_dist--;
	curNode = *sptr_node; 	  
	curDist = *sptr_dist;

	prefetch<PFHINT_L1>((mic_f*)node + 0);           
	prefetch<PFHINT_L1>((mic_f*)node + 1); 
	prefetch<PFHINT_L1>((mic_f*)node + 2); 
	prefetch<PFHINT_L1>((mic_f*)node + 3); 

#pragma unroll(4)
	for (unsigned int i=0; i<4; i++)
	  {
	    const BVH4i::NodeRef child = node->lower[i].child;

	    const mic_f lower_x =  one_time * nodeMB->lower[i].x + time * nodeMB->lower_t1[i].x;
	    const mic_f lower_y =  one_time * nodeMB->lower[i].y + time * nodeMB->lower_t1[i].y;
	    const mic_f lower_z =  one_time * nodeMB->lower[i].z + time * nodeMB->lower_t1[i].z;
	    const mic_f upper_x =  one_time * nodeMB->upper[i].x + time * nodeMB->upper_t1[i].x;
	    const mic_f upper_y =  one_time * nodeMB->upper[i].y + time * nodeMB->upper_t1[i].y;
	    const mic_f upper_z =  one_time * nodeMB->upper[i].z + time * nodeMB->upper_t1[i].z;

	    if (unlikely(i >=2 && child == BVH4i::invalidNode)) break;

	    const mic_f lclipMinX = msub(lower_x,rdir.x,org_rdir.x);
	    const mic_f lclipMinY = msub(lower_y,rdir.y,org_rdir.y);
	    const mic_f lclipMinZ = msub(lower_z,rdir.z,org_rdir.z);
	    const mic_f lclipMaxX = msub(upper_x,rdir.x,org_rdir.x);
	    const mic_f lclipMaxY = msub(upper_y,rdir.y,org_rdir.y);
	    const mic_f lclipMaxZ = msub(upper_z,rdir.z,org_rdir.z);
	    
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

  __forceinline void traverse_chunk_occluded(BVH4i::NodeRef &curNode,
					     mic_f &curDist,
					     const mic3f &rdir,
					     const mic3f &org_rdir,
					     const mic_f &ray_tnear,
					     const mic_f &ray_tfar,
					     const mic_m &m_active,
					     const mic_f &time,
					     BVH4i::NodeRef *__restrict__ &sptr_node,
					     mic_f *__restrict__ &sptr_dist,
					     const BVH4i::Node      * __restrict__ const nodes,
					     const unsigned int leaf_mask)
  {
    const mic_f one_time = (mic_f::one() - time);
    while (1)
      {
	/* test if this is a leaf node */
	if (unlikely(curNode.isLeaf(leaf_mask))) break;
          
	STAT3(shadow.trav_nodes,1,popcnt(ray_tfar > curDist),16);
	const BVH4i::Node* __restrict__ const node = curNode.node(nodes);
	const BVH4mb::Node* __restrict__ const nodeMB = (BVH4mb::Node*)node;
          
	prefetch<PFHINT_L1>((mic_f*)node + 0); 
	prefetch<PFHINT_L1>((mic_f*)node + 1); 
	prefetch<PFHINT_L1>((mic_f*)node + 2); 
	prefetch<PFHINT_L1>((mic_f*)node + 3); 

	/* pop of next node */
	sptr_node--;
	sptr_dist--;
	curNode = *sptr_node; 
	curDist = *sptr_dist;
          	 
#pragma unroll(4)
	for (unsigned int i=0; i<4; i++)
          {
	    const BVH4i::NodeRef child = node->lower[i].child;

	    const mic_f lower_x =  one_time * nodeMB->lower[i].x + time * nodeMB->lower_t1[i].x;
	    const mic_f lower_y =  one_time * nodeMB->lower[i].y + time * nodeMB->lower_t1[i].y;
	    const mic_f lower_z =  one_time * nodeMB->lower[i].z + time * nodeMB->lower_t1[i].z;
	    const mic_f upper_x =  one_time * nodeMB->upper[i].x + time * nodeMB->upper_t1[i].x;
	    const mic_f upper_y =  one_time * nodeMB->upper[i].y + time * nodeMB->upper_t1[i].y;
	    const mic_f upper_z =  one_time * nodeMB->upper[i].z + time * nodeMB->upper_t1[i].z;

	    if (unlikely(i >=2 && child == BVH4i::invalidNode)) break;

            const mic_f lclipMinX = msub(lower_x,rdir.x,org_rdir.x);
            const mic_f lclipMinY = msub(lower_y,rdir.y,org_rdir.y);
            const mic_f lclipMinZ = msub(lower_z,rdir.z,org_rdir.z);
            const mic_f lclipMaxX = msub(upper_x,rdir.x,org_rdir.x);
            const mic_f lclipMaxY = msub(upper_y,rdir.y,org_rdir.y);
            const mic_f lclipMaxZ = msub(upper_z,rdir.z,org_rdir.z);

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

};
