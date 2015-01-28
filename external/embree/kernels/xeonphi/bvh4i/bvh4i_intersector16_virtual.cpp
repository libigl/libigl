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

#include "bvh4i_intersector16_virtual.h"
#include "geometry/virtual_accel_intersector16.h"
#include "geometry/triangle1.h"
#include "geometry/virtual_accel_intersector1.h"

namespace embree
{
  namespace isa
  {
    static unsigned int BVH4I_LEAF_MASK = BVH4i::leaf_mask; // needed due to compiler efficiency bug

    template<typename VirtualGeometryIntersector16>
    void BVH4iIntersector16Virtual<VirtualGeometryIntersector16>::intersect(mic_i* valid_i, BVH4i* bvh, Ray16& ray)
    {
      /* near and node stack */
      __aligned(64) mic_f   stack_dist[3*BVH4i::maxDepth+1];
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

      if (unlikely(bvh == NULL)) return;

      /* load ray */
      const mic_m valid0   = *(mic_i*)valid_i != mic_i(0);
      const mic3f rdir     = rcp_safe(ray.dir);
 
      const mic3f org_rdir = ray.org * rdir;
      mic_f ray_tnear      = select(valid0,ray.tnear,pos_inf);
      mic_f ray_tfar       = select(valid0,ray.tfar ,neg_inf);
      const mic_f inf      = mic_f(pos_inf);
      
      /* allocate stack and push root node */
      stack_node[0] = BVH4i::invalidNode;
      stack_dist[0] = inf;
      stack_node[1] = bvh->root;
      stack_dist[1] = ray_tnear; 
      NodeRef* __restrict__ sptr_node = stack_node + 2;
      mic_f*   __restrict__ sptr_dist = stack_dist + 2;
      
      const Node     * __restrict__ nodes  = (Node    *)bvh->nodePtr();
      const Triangle1 * __restrict__ accel = (Triangle1*)bvh->triPtr();

      if (unlikely(accel == NULL)) return;

      while (1)
      {
        /* pop next node from stack */
        NodeRef curNode = *(sptr_node-1);
        mic_f curDist   = *(sptr_dist-1);
        sptr_node--;
        sptr_dist--;
	const mic_m m_stackDist = ray_tfar > curDist;

	/* stack emppty ? */
        if (unlikely(curNode == BVH4i::invalidNode))  break;
        
        /* cull node if behind closest hit point */
        if (unlikely(none(m_stackDist))) continue;
	        
	const unsigned int leaf_mask = BVH4I_LEAF_MASK;

        while (1)
        {
          /* test if this is a leaf node */
          if (unlikely(curNode.isLeaf(leaf_mask))) break;
          
          STAT3(normal.trav_nodes,1,popcnt(ray_tfar > curDist),16);
          const Node* __restrict__ const node = curNode.node(nodes);

          /* pop of next node */
          sptr_node--;
          sptr_dist--;
          curNode = *sptr_node; 
          curDist = *sptr_dist;
          
	  prefetch<PFHINT_L1>((mic_f*)node + 1); // depth first order, prefetch		

#pragma unroll(4)
          for (unsigned int i=0; i<4; i++)
          {
            //const NodeRef child = node->children[i];
	    const NodeRef child = node->lower[i].child;

            //if (unlikely(child == BVH4i::emptyNode)) break;
	    
            const mic_f lclipMinX = msub(node->lower[i].x,rdir.x,org_rdir.x);
            const mic_f lclipMinY = msub(node->lower[i].y,rdir.y,org_rdir.y);
            const mic_f lclipMinZ = msub(node->lower[i].z,rdir.z,org_rdir.z);
            const mic_f lclipMaxX = msub(node->upper[i].x,rdir.x,org_rdir.x);
            const mic_f lclipMaxY = msub(node->upper[i].y,rdir.y,org_rdir.y);
            const mic_f lclipMaxZ = msub(node->upper[i].z,rdir.z,org_rdir.z);
	    
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
              assert(sptr_node - stack_node < BVH4i::maxDepth);
            }	      
          }
        }
        
        /* return if stack is empty */
        if (unlikely(curNode == BVH4i::invalidNode)) break;
        
        /* intersect leaf */
        const mic_m valid_leaf = ray_tfar > curDist;
        STAT3(normal.trav_leaves,1,popcnt(valid_leaf),16);
 
	//////////////////////////////////////////////////////////////////////////////////////////////////

	unsigned int items = curNode.items();
	unsigned int index = (curNode.offset() >> 6); /* array of AccelSetItems */
	Primitive *accel_ptr = (Primitive*)accel + index;

        VirtualGeometryIntersector16::intersect(valid_leaf,ray,accel_ptr,items,bvh->geometry);
	//////////////////////////////////////////////////////////////////////////////////////////////////

        ray_tfar = select(valid_leaf,ray.tfar,ray_tfar);
      }
    }
    
    template<typename VirtualGeometryIntersector16>
    void BVH4iIntersector16Virtual<VirtualGeometryIntersector16>::occluded(mic_i* valid_i, BVH4i* bvh, Ray16& ray)
    {
      /* allocate stack */
      __aligned(64) mic_f    stack_dist[3*BVH4i::maxDepth+1];
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

      if (unlikely(bvh == NULL)) return;

      /* load ray */
      const mic_m valid = *(mic_i*)valid_i != mic_i(0);
      mic_m m_terminated = !valid;
      const mic3f rdir = rcp_safe(ray.dir);
      const mic3f org_rdir = ray.org * rdir;
      mic_f ray_tnear = select(valid,ray.tnear,pos_inf);
      mic_f ray_tfar  = select(valid,ray.tfar ,neg_inf);
      const mic_f inf = mic_f(pos_inf);
      
      /* push root node */
      stack_node[0] = BVH4i::invalidNode;
      stack_dist[0] = inf;
      stack_node[1] = bvh->root;
      stack_dist[1] = ray_tnear; 
      NodeRef* __restrict__ sptr_node = stack_node + 2;
      mic_f*   __restrict__ sptr_dist = stack_dist + 2;
      
      const Node     * __restrict__ nodes = (Node    *)bvh->nodePtr();
      const Triangle1 * __restrict__ accel = (Triangle1*)bvh->triPtr();

      if (unlikely(accel == NULL)) return;

      while (1)
      {
	const mic_m m_active = !m_terminated;

        /* pop next node from stack */
        NodeRef curNode = *(sptr_node-1);
        mic_f curDist   = *(sptr_dist-1);
        sptr_node--;
        sptr_dist--;
	const mic_m m_stackDist = gt(m_active,ray_tfar,curDist);

	/* stack emppty ? */
        if (unlikely(curNode == BVH4i::invalidNode))  break;
        
        /* cull node if behind closest hit point */
        if (unlikely(none(m_stackDist))) continue;
	
        while (1)
        {
          /* test if this is a leaf node */
          if (unlikely(curNode.isLeaf())) break;
          
          STAT3(shadow.trav_nodes,1,popcnt(ray_tfar > curDist),16);
          const Node* __restrict__ const node = curNode.node(nodes);
          
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
	    const NodeRef child = node->lower[i].child;

            //if (unlikely(child == BVH4i::emptyNode)) break;
            
            const mic_f lclipMinX = msub(node->lower[i].x,rdir.x,org_rdir.x);
            const mic_f lclipMinY = msub(node->lower[i].y,rdir.y,org_rdir.y);
            const mic_f lclipMinZ = msub(node->lower[i].z,rdir.z,org_rdir.z);
            const mic_f lclipMaxX = msub(node->upper[i].x,rdir.x,org_rdir.x);
            const mic_f lclipMaxY = msub(node->upper[i].y,rdir.y,org_rdir.y);
            const mic_f lclipMaxZ = msub(node->upper[i].z,rdir.z,org_rdir.z);	    

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
              assert(sptr_node - stack_node < BVH4i::maxDepth);
            }	      
          }
        }
        
        /* return if stack is empty */
        if (unlikely(curNode == BVH4i::invalidNode)) break;
        
        /* intersect leaf */
        mic_m valid_leaf = gt(m_active,ray_tfar,curDist);
        STAT3(shadow.trav_leaves,1,popcnt(valid_leaf),16);

	//////////////////////////////////////////////////////////////////////////////////////////////////

	unsigned int items = curNode.items();
	unsigned int index = (curNode.offset() >> 6); /* array of AccelSetItems */
	Primitive *accel_ptr = (Primitive *)accel + index;

        m_terminated |= valid_leaf & VirtualGeometryIntersector16::occluded(valid_leaf,ray,accel_ptr,items,bvh->geometry);

	//////////////////////////////////////////////////////////////////////////////////////////////////

        if (unlikely(all(m_terminated))) break;
        ray_tfar = select(m_terminated,neg_inf,ray_tfar);
      }
      store16i(valid & m_terminated,&ray.geomID,0);
    }

    template<typename VirtualGeometryIntersector1>
    void BVH4iIntersector1Virtual<VirtualGeometryIntersector1>::intersect(BVH4i* bvh, Ray& ray)
    {
      /* near and node stack */
      __aligned(64) float   stack_dist[3*BVH4i::maxDepth+1];
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

      if (unlikely(bvh == NULL)) return;

      /* setup */
      const mic3f rdir16     = rcp_safe(mic3f(mic_f(ray.dir.x),mic_f(ray.dir.y),mic_f(ray.dir.z)));
      const mic_f inf        = mic_f(pos_inf);
      const mic_f zero       = mic_f::zero();

      store16f(stack_dist,inf);

      const Node      * __restrict__ nodes = (Node    *)bvh->nodePtr();
      const Triangle1 * __restrict__ accel = (Triangle1*)bvh->triPtr();

      if (unlikely(accel == NULL)) return;

      stack_node[0] = BVH4i::invalidNode;      
      stack_node[1] = bvh->root;

      size_t sindex = 2;

      const mic_f org_xyz      = loadAOS4to16f(ray.org.x,ray.org.y,ray.org.z);
      const mic_f dir_xyz      = loadAOS4to16f(ray.dir.x,ray.dir.y,ray.dir.z);
      const mic_f rdir_xyz     = loadAOS4to16f(rdir16.x[0],rdir16.y[0],rdir16.z[0]);
      const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
      const mic_f min_dist_xyz = broadcast1to16f(&ray.tnear);
      mic_f       max_dist_xyz = broadcast1to16f(&ray.tfar);
	  
      const unsigned int leaf_mask = BVH4I_LEAF_MASK;
	  
      while (1)
	{
	  NodeRef curNode = stack_node[sindex-1];
	  sindex--;
            
	  while (1) 
	    {
	      /* test if this is a leaf node */
	      if (unlikely(curNode.isLeaf(leaf_mask))) break;
        
	      const Node* __restrict__ const node = curNode.node(nodes);	      
	      const float* __restrict const plower = (float*)node->lower;
	      const float* __restrict const pupper = (float*)node->upper;

	      
	      prefetch<PFHINT_L1>((char*)node + 0);
	      prefetch<PFHINT_L1>((char*)node + 64);
        
	      /* intersect single ray with 4 bounding boxes */
	      const mic_f tLowerXYZ = load16f(plower) * rdir_xyz - org_rdir_xyz;
	      const mic_f tUpperXYZ = load16f(pupper) * rdir_xyz - org_rdir_xyz;
	      const mic_f tLower = mask_min(0x7777,min_dist_xyz,tLowerXYZ,tUpperXYZ);
	      const mic_f tUpper = mask_max(0x7777,max_dist_xyz,tLowerXYZ,tUpperXYZ);

	      sindex--;
	      curNode = stack_node[sindex]; // early pop of next node

	      const Node* __restrict__ const next = curNode.node(nodes);
	      prefetch<PFHINT_L2>((char*)next + 0);
	      prefetch<PFHINT_L2>((char*)next + 64);

	      const mic_f tNear = vreduce_max4(tLower);
	      const mic_f tFar  = vreduce_min4(tUpper);  
	      const mic_m hitm = le(0x8888,tNear,tFar);
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
        
	      const mic_m closest_child = eq(hitm,min_dist,tNear);
	      const unsigned long closest_child_pos = bitscan64(closest_child);
	      const mic_m m_pos = andn(hitm,andn(closest_child,(mic_m)((unsigned int)closest_child - 1)));
	      const mic_i plower_node = load16i((int*)plower);
	      compactustore16i(m_pos,&stack_node[old_sindex],plower_node);
	      curNode = ((unsigned int*)plower)[closest_child_pos];
	      compactustore16f(m_pos,&stack_dist[old_sindex],tNear);
	    }
	  	    

	  /* return if stack is empty */
	  if (unlikely(curNode == BVH4i::invalidNode)) break;

	  //////////////////////////////////////////////////////////////////////////////////////////////////

	  unsigned int items = curNode.items();
	  unsigned int index = (curNode.offset() >> 6); /* array of AccelSetItems */
	  Primitive *accel_ptr = (Primitive *)accel + index;

	  VirtualGeometryIntersector1::intersect(ray,accel_ptr,items,bvh->geometry);

	  if (unlikely(any(max_dist_xyz != broadcast1to16f(&ray.tfar))))
	    {
	      max_dist_xyz = broadcast1to16f(&ray.tfar);

	      /* compact the stack if size of stack >= 2 */
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
		      sindex = countbits(m_stack_compact_low) + countbits(m_stack_compact_high);
		      assert ((unsigned int)m_num_stack_high == ((mic_m::shift1[sindex] - 1) >> 16));
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
		} // sindex
	    }

	  //////////////////////////////////////////////////////////////////////////////////////////////////

	}	  
    }

    template<typename VirtualGeometryIntersector1>
    void BVH4iIntersector1Virtual<VirtualGeometryIntersector1>::occluded(BVH4i* bvh, Ray& ray)
    {
      /* near and node stack */
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

      if (unlikely(bvh == NULL)) return;

      /* setup */
      const mic3f rdir16      = rcp_safe(mic3f(ray.dir.x,ray.dir.y,ray.dir.z));
      const mic_f inf         = mic_f(pos_inf);
      const mic_f zero        = mic_f::zero();

      const Node      * __restrict__ nodes = (Node     *)bvh->nodePtr();
      const Triangle1 * __restrict__ accel = (Triangle1*)bvh->triPtr();

      if (unlikely(accel == NULL)) return;

      stack_node[0] = BVH4i::invalidNode;
      stack_node[1] = bvh->root;
      size_t sindex = 2;

      const mic_f org_xyz      = loadAOS4to16f(ray.org.x,ray.org.y,ray.org.z);
      const mic_f dir_xyz      = loadAOS4to16f(ray.dir.x,ray.dir.y,ray.dir.z);
      const mic_f rdir_xyz     = loadAOS4to16f(rdir16.x[0],rdir16.y[0],rdir16.z[0]);
      const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
      const mic_f min_dist_xyz = broadcast1to16f(&ray.tnear);
      const mic_f max_dist_xyz = broadcast1to16f(&ray.tfar);

      const unsigned int leaf_mask = BVH4I_LEAF_MASK;
	  
      while (1)
	{
	  NodeRef curNode = stack_node[sindex-1];
	  sindex--;
            
	  while (1) 
	    {
	      /* test if this is a leaf node */
	      if (unlikely(curNode.isLeaf(leaf_mask))) break;
        
	      const Node* __restrict__ const node = curNode.node(nodes);
	      const float* __restrict const plower = (float*)node->lower;
	      const float* __restrict const pupper = (float*)node->upper;

	      prefetch<PFHINT_L1>((char*)node + 0);
	      prefetch<PFHINT_L1>((char*)node + 64);
        
	      /* intersect single ray with 4 bounding boxes */
	      const mic_f tLowerXYZ = load16f(plower) * rdir_xyz - org_rdir_xyz;
	      const mic_f tUpperXYZ = load16f(pupper) * rdir_xyz - org_rdir_xyz;
	      const mic_f tLower = mask_min(0x7777,min_dist_xyz,tLowerXYZ,tUpperXYZ);
	      const mic_f tUpper = mask_max(0x7777,max_dist_xyz,tLowerXYZ,tUpperXYZ);

	      sindex--;
	      curNode = stack_node[sindex]; 

	      const Node* __restrict__ const next = curNode.node(nodes);
	      prefetch<PFHINT_L2>((char*)next + 0);
	      prefetch<PFHINT_L2>((char*)next + 64);

	      const mic_f tNear = vreduce_max4(tLower);
	      const mic_f tFar  = vreduce_min4(tUpper);  
	      const mic_m hitm = le(0x8888,tNear,tFar);
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
	      const unsigned old_sindex = sindex;
	      sindex += countbits(hiti) - 1;
	      assert(sindex < 3*BVH4i::maxDepth+1);
        
	      const mic_m closest_child = eq(hitm,min_dist,tNear);
	      const unsigned long closest_child_pos = bitscan64(closest_child);
	      const mic_m m_pos = andn(hitm,andn(closest_child,(mic_m)((unsigned int)closest_child - 1)));
	      const mic_i plower_node = load16i((int*)plower);
	      curNode = ((unsigned int*)plower)[closest_child_pos];
	      compactustore16i(m_pos,&stack_node[old_sindex],plower_node);
	    }
	  
	    

	  /* return if stack is empty */
	  if (unlikely(curNode == BVH4i::invalidNode)) break;



	  //////////////////////////////////////////////////////////////////////////////////////////////////

	  unsigned int items = curNode.items();
	  unsigned int index = (curNode.offset() >> 6); /* array of AccelSetItems */
	  Primitive *accel_ptr = (Primitive *)accel + index;

	  if (VirtualGeometryIntersector1::occluded(ray,accel_ptr,items,bvh->geometry)) {
	    ray.geomID = 0;
	    return;
	  }

	  //////////////////////////////////////////////////////////////////////////////////////////////////

	}
    }
    
    
    DEFINE_INTERSECTOR16   (BVH4iVirtualGeometryIntersector16, BVH4iIntersector16Virtual<VirtualAccelIntersector16>);
    DEFINE_INTERSECTOR1    (BVH4iVirtualGeometryIntersector1, BVH4iIntersector1Virtual<VirtualAccelIntersector1>);

  }
}
