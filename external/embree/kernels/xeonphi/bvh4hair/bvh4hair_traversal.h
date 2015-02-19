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

#include "bvh4hair.h"

namespace embree
{
  //TODO: pre-multiply v with 1/127.0f

  static __forceinline mic_f xfm_row_vector(const mic_f &row0,
					    const mic_f &row1,
					    const mic_f &row2,
					    const mic_f &v)
  {
    const mic_f x = swAAAA(v);
    const mic_f y = swBBBB(v);
    const mic_f z = swCCCC(v);
    const mic_f ret = x * row0 + y * row1 + z * row2;
    return ret;
  }

  static __forceinline mic_f xfm_row_point(const mic_f &row0,
					   const mic_f &row1,
					   const mic_f &row2,
					   const mic_f &v)
  {
    const mic_f x = swAAAA(v);
    const mic_f y = swBBBB(v);
    const mic_f z = swCCCC(v);
    const mic_f trans  = select(0x4444,swDDDD(row2),select(0x2222,swDDDD(row1),swDDDD(row0)));   
    const mic_f ret = x * row0 + y * row1 + z * row2 + trans;
    return ret;
  }

  static __forceinline mic_f xfm(const mic_f &v, const BVH4Hair::UnalignedNode &node)
  {
    // alternative: 'rows' at 'columns' for lane dot products

    const mic_f x = ldot3_xyz(v,node.getRow(0));
    const mic_f y = ldot3_xyz(v,node.getRow(1));
    const mic_f z = ldot3_xyz(v,node.getRow(2));
    const mic_f ret = select(0x4444,z,select(0x2222,y,x));
    return ret;
  }

  __forceinline void traverse_single_intersect(BVH4Hair::NodeRef &curNode,
					       size_t &sindex,
					       const mic_f &dir_xyz,
					       const mic_f &org_xyz1,
					       const mic_f &rdir_xyz,
					       const mic_f &org_rdir_xyz,
					       const mic_f &min_dist_xyz,
					       const mic_f &max_dist_xyz,
					       BVH4Hair::NodeRef *__restrict__ const stack_node,
					       float   *__restrict__ const stack_dist,
					       const void      * __restrict__ const nodes,
					       const unsigned int leaf_mask,
					       const unsigned int alignednode_mask)
  {
    const mic_m m7777 = 0x7777; 

    const mic_i invalidNode = mic_i::neg_one();

    while (1) 
      {
	if (unlikely(curNode.isLeaf(leaf_mask))) break;

	mic_m hitm;
	mic_f tNear,tFar;

	const BVH4Hair::UnalignedNode *__restrict__ const u_node = (BVH4Hair::UnalignedNode *)curNode.node(nodes);

	BVH4Hair::NodeRef *__restrict__ ref; 
	
	STAT3(normal.trav_nodes,1,1,1);	    

	if (likely((curNode & alignednode_mask) == 0)) // unaligned nodes
	  {	    
	    u_node->prefetchNode<PFHINT_L1>();
	    
	    const mic_f row0 = u_node->getRow(0);
	    const mic_f row1 = u_node->getRow(1);
	    const mic_f row2 = u_node->getRow(2);
		  	
	    const mic_f xfm_dir_xyz = xfm_row_vector(row0,row1,row2,dir_xyz);
	    const mic_f xfm_org_xyz = xfm_row_vector(row0,row1,row2,org_xyz1); 

	    const mic_f rcp_xfm_dir_xyz = rcp_safe( xfm_dir_xyz );
		
	    const mic_f xfm_org_rdir_xyz = xfm_org_xyz * rcp_xfm_dir_xyz;
	    const mic_f tLowerXYZ = msub(rcp_xfm_dir_xyz,load16f(u_node->lower),xfm_org_rdir_xyz);
	    const mic_f tUpperXYZ = msub(rcp_xfm_dir_xyz,load16f(u_node->upper),xfm_org_rdir_xyz);
		   
	    const mic_f tLower = select(m7777,min(tLowerXYZ,tUpperXYZ),min_dist_xyz);
	    const mic_f tUpper = select(m7777,max(tLowerXYZ,tUpperXYZ),max_dist_xyz);
	    ref = u_node->nodeRefPtr();

	    hitm = ne(0x8888,invalidNode,load16i((const int*)ref));

	    /* early pop of next node */
	    sindex--;
	    curNode = stack_node[sindex];

	    tNear = vreduce_max4(tLower);
	    tFar  = vreduce_min4(tUpper);  


	  }
	else
	  {
	    prefetch<PFHINT_L1>((mic_f*)u_node + 0);
	    prefetch<PFHINT_L1>((mic_f*)u_node + 1);

	    const BVH4Hair::AlignedNode* __restrict__ node = (BVH4Hair::AlignedNode*)((size_t)u_node ^ BVH4Hair::alignednode_mask);

	    //node->prefetchNode<PFHINT_L1>();

	    const float* __restrict const plower = (float*)node->lower;
	    const float* __restrict const pupper = (float*)node->upper;

        
	    /* intersect single ray with 4 bounding boxes */
	    ref = node->nodeRefPtr();
	    hitm = ne(0x8888,invalidNode,load16i((const int*)ref));


	    const mic_f tLowerXYZ = load16f(plower) * rdir_xyz - org_rdir_xyz;
	    const mic_f tUpperXYZ = load16f(pupper) * rdir_xyz - org_rdir_xyz;
	    const mic_f tLower = mask_min(0x7777,min_dist_xyz,tLowerXYZ,tUpperXYZ);
	    const mic_f tUpper = mask_max(0x7777,max_dist_xyz,tLowerXYZ,tUpperXYZ);


	    /* if no child is hit, continue with early popped child */	    
	    sindex--;
	    curNode = stack_node[sindex];
	    tNear = vreduce_max4(tLower);
	    tFar  = vreduce_min4(tUpper);  
	  }

	hitm = le(hitm,tNear,tFar);
	const mic_f tNear_pos = select(hitm,tNear,inf);

	STAT3(normal.trav_hit_boxes[countbits(hitm)],1,1,1);

	/* if no child is hit, continue with early popped child */
	if (unlikely(none(hitm))) continue;

		  
	sindex++;        
	const unsigned long hiti = toInt(hitm);
	const unsigned long pos_first = bitscan64(hiti);
	const unsigned long num_hitm = countbits(hiti); 
        
	assert(num_hitm <= 4);

	/* if a single child is hit, continue with that child */
	//curNode = u_node->child_ref(pos_first);
	curNode = ref[pos_first];

	assert(curNode != BVH4Hair::invalidNode);

	if (likely(num_hitm == 1)) continue;
        

	/* if two children are hit, push in correct order */
	const unsigned long pos_second = bitscan64(pos_first,hiti);
	if (likely(num_hitm == 2))
	  {
	    const unsigned int dist_first  = ((unsigned int*)&tNear)[pos_first];
	    const unsigned int dist_second = ((unsigned int*)&tNear)[pos_second];
	    const BVH4Hair::NodeRef node_first  = curNode;
	    const BVH4Hair::NodeRef node_second = ref[pos_second];

	    assert(node_first  != BVH4Hair::invalidNode);
	    assert(node_second != BVH4Hair::invalidNode);
          
	    if (dist_first <= dist_second)
	      {			  
		stack_node[sindex] = node_second;
		((unsigned int*)stack_dist)[sindex] = dist_second;                      
		sindex++;
		assert(sindex < 3*BVH4Hair::maxDepth+1);
		continue;
	      }
	    else
	      {
		stack_node[sindex] = node_first;
		((unsigned int*)stack_dist)[sindex] = dist_first;
		curNode = node_second;
		sindex++;
		assert(sindex < 3*BVH4Hair::maxDepth+1);
		continue;
	      }
	  }

	/* continue with closest child and push all others */

	const mic_f min_dist = set_min_lanes(tNear_pos);
	const unsigned int old_sindex = sindex;
	sindex += countbits(hiti) - 1;
	assert(sindex < 3*BVH4Hair::maxDepth+1);
	const mic_i children = load16i((const int*)ref); 
        
	const mic_m closest_child = eq(hitm,min_dist,tNear);
	const unsigned long closest_child_pos = bitscan64(closest_child);
	const mic_m m_pos = andn(hitm,andn(closest_child,(mic_m)((unsigned int)closest_child - 1)));
	curNode = ref[closest_child_pos];

	compactustore16f(m_pos,&stack_dist[old_sindex],tNear);
	compactustore16i(m_pos,&stack_node[old_sindex],children);

	if (unlikely(((unsigned int*)stack_dist)[sindex-3] < ((unsigned int*)stack_dist)[sindex-2]))
	  {
	    std::swap(((unsigned int*)stack_dist)[sindex-2],((unsigned int*)stack_dist)[sindex-3]);
	    std::swap(((unsigned int*)stack_node)[sindex-2],((unsigned int*)stack_node)[sindex-3]);
	  }

	if (unlikely(((unsigned int*)stack_dist)[sindex-2] < ((unsigned int*)stack_dist)[sindex-1]))
	  {
	    std::swap(((unsigned int*)stack_dist)[sindex-1],((unsigned int*)stack_dist)[sindex-2]);
	    std::swap(((unsigned int*)stack_node)[sindex-1],((unsigned int*)stack_node)[sindex-2]);
	  }
      }
 
  }



  __forceinline void traverse_single_occluded(BVH4Hair::NodeRef &curNode,
					      size_t &sindex,
					      const mic_f &dir_xyz,
					      const mic_f &org_xyz1,
					      const mic_f &rdir_xyz,
					      const mic_f &org_rdir_xyz,
					      const mic_f &min_dist_xyz,
					      const mic_f &max_dist_xyz,
					      BVH4Hair::NodeRef *__restrict__ const stack_node,
					      const void      * __restrict__ const nodes,
					      const unsigned int leaf_mask)
  {
    const mic_m m7777 = 0x7777; 
    const mic_i invalidNode = mic_i::neg_one();
    
    while (1) 
      {
	if (unlikely(curNode.isLeaf(leaf_mask))) break;
		  
	mic_m hitm;
	mic_f tNear,tFar;

	const BVH4Hair::UnalignedNode *__restrict__ const u_node = (BVH4Hair::UnalignedNode *)curNode.node(nodes);

	BVH4Hair::NodeRef *__restrict__ ref; 
	
	STAT3(normal.trav_nodes,1,1,1);	    

	if (likely(((size_t)curNode & BVH4Hair::alignednode_mask) == 0)) // unaligned nodes
	  {	    
	    u_node->prefetchNode<PFHINT_L1>();
	    
	    const mic_f row0 = u_node->getRow(0);
	    const mic_f row1 = u_node->getRow(1);
	    const mic_f row2 = u_node->getRow(2);
		  	
	    const mic_f xfm_dir_xyz = xfm_row_vector(row0,row1,row2,dir_xyz);
	    const mic_f xfm_org_xyz = xfm_row_vector(row0,row1,row2,org_xyz1); 

	    const mic_f rcp_xfm_dir_xyz = rcp_safe( xfm_dir_xyz );
		
	    const mic_f xfm_org_rdir_xyz = xfm_org_xyz * rcp_xfm_dir_xyz;
	    const mic_f tLowerXYZ = msub(rcp_xfm_dir_xyz,load16f(u_node->lower),xfm_org_rdir_xyz);
	    const mic_f tUpperXYZ = msub(rcp_xfm_dir_xyz,load16f(u_node->upper),xfm_org_rdir_xyz);
		   
	    const mic_f tLower = select(m7777,min(tLowerXYZ,tUpperXYZ),min_dist_xyz);
	    const mic_f tUpper = select(m7777,max(tLowerXYZ,tUpperXYZ),max_dist_xyz);

	    /* early pop of next node */
	    sindex--;
	    curNode = stack_node[sindex];

	    tNear = vreduce_max4(tLower);
	    tFar  = vreduce_min4(tUpper);  
	    ref = u_node->nodeRefPtr();

	    hitm = ne(0x8888,invalidNode,load16i((const int*)ref));

	  }
	else
	  {
	    prefetch<PFHINT_L1>((mic_f*)u_node + 0);
	    prefetch<PFHINT_L1>((mic_f*)u_node + 1);

	    const BVH4Hair::AlignedNode* __restrict__ node = (BVH4Hair::AlignedNode*)((size_t)u_node ^ BVH4Hair::alignednode_mask);

	    //node->prefetchNode<PFHINT_L1>();

	    const float* __restrict const plower = (float*)node->lower;
	    const float* __restrict const pupper = (float*)node->upper;

        
	    /* intersect single ray with 4 bounding boxes */
	    ref = node->nodeRefPtr();
	    hitm = ne(0x8888,invalidNode,load16i((const int*)ref));


	    const mic_f tLowerXYZ = load16f(plower) * rdir_xyz - org_rdir_xyz;
	    const mic_f tUpperXYZ = load16f(pupper) * rdir_xyz - org_rdir_xyz;
	    const mic_f tLower = mask_min(0x7777,min_dist_xyz,tLowerXYZ,tUpperXYZ);
	    const mic_f tUpper = mask_max(0x7777,max_dist_xyz,tLowerXYZ,tUpperXYZ);


	    /* if no child is hit, continue with early popped child */	    
	    sindex--;
	    curNode = stack_node[sindex];
	    tNear = vreduce_max4(tLower);
	    tFar  = vreduce_min4(tUpper);  
	  }

	hitm = le(hitm,tNear,tFar);
	const mic_f tNear_pos = select(hitm,tNear,inf);

	STAT3(normal.trav_hit_boxes[countbits(hitm)],1,1,1);

	/* if no child is hit, continue with early popped child */
	if (unlikely(none(hitm))) continue;

		  
	sindex++;        
	const unsigned long hiti = toInt(hitm);
	const unsigned long pos_first = bitscan64(hiti);
	const unsigned long num_hitm = countbits(hiti); 
        
	assert(num_hitm <= 4);

	/* if a single child is hit, continue with that child */
	//curNode = u_node->child_ref(pos_first);
	curNode = ref[pos_first];

	assert(curNode != BVH4Hair::invalidNode);

	if (likely(num_hitm == 1)) continue;
        

	/* if two children are hit, push in correct order */
	const unsigned long pos_second = bitscan64(pos_first,hiti);
	if (likely(num_hitm == 2))
	  {
	    const unsigned int dist_first  = ((unsigned int*)&tNear)[pos_first];
	    const unsigned int dist_second = ((unsigned int*)&tNear)[pos_second];
	    const BVH4Hair::NodeRef node_first  = curNode;
	    const BVH4Hair::NodeRef node_second = ref[pos_second];

	    assert(node_first  != BVH4Hair::invalidNode);
	    assert(node_second != BVH4Hair::invalidNode);
          
	    if (dist_first <= dist_second)
	      {			  
		stack_node[sindex] = node_second;
		sindex++;
		assert(sindex < 3*BVH4Hair::maxDepth+1);
		continue;
	      }
	    else
	      {
		stack_node[sindex] = node_first;
		curNode = node_second;
		sindex++;
		assert(sindex < 3*BVH4Hair::maxDepth+1);
		continue;
	      }
	  }

	/* continue with closest child and push all others */

	const mic_f min_dist = set_min_lanes(tNear_pos);
	const unsigned int old_sindex = sindex;
	sindex += countbits(hiti) - 1;
	assert(sindex < 3*BVH4Hair::maxDepth+1);
	const mic_i children = load16i((const int*)ref); 
        
	const mic_m closest_child = eq(hitm,min_dist,tNear);
	const unsigned long closest_child_pos = bitscan64(closest_child);
	const mic_m m_pos = andn(hitm,andn(closest_child,(mic_m)((unsigned int)closest_child - 1)));
	curNode = ref[closest_child_pos];

	compactustore16i(m_pos,&stack_node[old_sindex],children);
      }
  }



__forceinline void compactStack(BVH4Hair::NodeRef *__restrict__ const stack_node,
				float   *__restrict__ const stack_dist,
				size_t &sindex,
				const unsigned int current_dist,
				const mic_f &max_dist_xyz)
{
#if 0
  size_t new_sindex = 1;
  for (size_t s=1;s<sindex;s++)
    if (((unsigned int*)stack_dist)[s] <= current_dist)
      {
	stack_dist[new_sindex] = stack_dist[s];
	stack_node[new_sindex] = stack_node[s];
	new_sindex++;
      }
  sindex = new_sindex;
#else
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

#endif   
}
  
};

/* const mic_f xfm_dir_xyz = xfm_row_vector(row0,row1,row2,dir_xyz); */
/* const mic_f xfm_org_xyz = xfm_row_point (row0,row1,row2,org_xyz1); // only use xfm_row_vector */

/* const mic_f rcp_xfm_dir_xyz = rcp_safe( xfm_dir_xyz ); */

/* const mic_f tLowerXYZ = -(xfm_org_xyz * rcp_xfm_dir_xyz); */
/* const mic_f tUpperXYZ = (mic_f::one()  - xfm_org_xyz) * rcp_xfm_dir_xyz; */

/* const mic_f min_dist = set_min_lanes(tNear_pos); */
/* assert(sindex < 3*BVH4Hair::maxDepth+1); */
        
/* const mic_m closest_child = eq(hitm,min_dist,tNear); */
/* const unsigned long closest_child_pos = bitscan64(closest_child); */
/* const mic_m m_pos = andn(hitm,andn(closest_child,(mic_m)((unsigned int)closest_child - 1))); */


/* curNode = u_node->child_ref(closest_child_pos); */


/* assert(curNode  != BVH4Hair::invalidNode); */

/* long i = -1; */
/* while((i = bitscan64(i,m_pos)) != BITSCAN_NO_BIT_SET_64)	     */
/*   { */
/* 	((unsigned int*)stack_dist)[sindex] = ((unsigned int*)&tNear)[i];		       */
/* 	stack_node[sindex] = u_node->child_ref(i); */
/* 	assert(stack_node[sindex]  != BVH4Hair::invalidNode); */
/* 	sindex++; */
/*   } */
