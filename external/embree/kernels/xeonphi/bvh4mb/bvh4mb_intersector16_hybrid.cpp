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

#include "bvh4mb_intersector16_hybrid.h"
#include "bvh4mb_traversal.h"
#include "bvh4mb_leaf_intersector.h"

#define SWITCH_ON_DOWN_TRAVERSAL 1

namespace embree
{
  namespace isa
  {
    static unsigned int BVH4I_LEAF_MASK = BVH4i::leaf_mask; // needed due to compiler efficiency bug

    template<typename LeafIntersector>
    void BVH4mbIntersector16Hybrid<LeafIntersector>::intersect(mic_i* valid_i, BVH4mb* bvh, Ray16& ray16)
    {
      /* near and node stack */
      __aligned(64) mic_f   stack_dist[3*BVH4i::maxDepth+1];
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];
      __aligned(64) NodeRef stack_node_single[3*BVH4i::maxDepth+1]; 

      /* load ray */
      const mic_m valid0     = *(mic_i*)valid_i != mic_i(0);
      const mic3f rdir16     = rcp_safe(ray16.dir);
      const mic3f org_rdir16 = ray16.org * rdir16;
      mic_f ray_tnear        = select(valid0,ray16.tnear,pos_inf);
      mic_f ray_tfar         = select(valid0,ray16.tfar ,neg_inf);
      const mic_f inf        = mic_f(pos_inf);
      
      /* allocate stack and push root node */
      stack_node[0] = BVH4i::invalidNode;
      stack_dist[0] = inf;
      stack_node[1] = bvh->root;
      stack_dist[1] = ray_tnear; 
      NodeRef* __restrict__ sptr_node = stack_node + 2;
      mic_f*   __restrict__ sptr_dist = stack_dist + 2;
      
      const Node      * __restrict__ nodes = (Node     *)bvh->nodePtr();
      const BVH4mb::Triangle01 * __restrict__ accel = (BVH4mb::Triangle01 *)bvh->triPtr();

      while (1) pop:
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
        
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/* switch to single ray mode */
        if (unlikely(countbits(m_stackDist) <= BVH4i::hybridSIMDUtilSwitchThreshold)) 
	  {
	    float   *__restrict__ stack_dist_single = (float*)sptr_dist;
	    store16f(stack_dist_single,inf);

	    /* traverse single ray */	  	  
	    long rayIndex = -1;
	    while((rayIndex = bitscan64(rayIndex,m_stackDist)) != BITSCAN_NO_BIT_SET_64) 
	      {	    
		stack_node_single[0] = BVH4i::invalidNode;
		stack_node_single[1] = curNode;
		size_t sindex = 2;

		const mic_f org_xyz      = loadAOS4to16f(rayIndex,ray16.org.x,ray16.org.y,ray16.org.z);
		const mic_f dir_xyz      = loadAOS4to16f(rayIndex,ray16.dir.x,ray16.dir.y,ray16.dir.z);
		const mic_f rdir_xyz     = loadAOS4to16f(rayIndex,rdir16.x,rdir16.y,rdir16.z);
		const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
		const mic_f min_dist_xyz = broadcast1to16f(&ray16.tnear[rayIndex]);
		mic_f       max_dist_xyz = broadcast1to16f(&ray16.tfar[rayIndex]);
		const mic_f time         = broadcast1to16f(&ray16.time[rayIndex]);

		const unsigned int leaf_mask = BVH4I_LEAF_MASK;

		while (1) 
		  {
		    NodeRef curNode = stack_node_single[sindex-1];
		    sindex--;
            
		    traverse_single_intersect(curNode,
					      sindex,
					      rdir_xyz,
					      org_rdir_xyz,
					      min_dist_xyz,
					      max_dist_xyz,
					      time,
					      stack_node_single,
					      stack_dist_single,
					      nodes,
					      leaf_mask);	    

		    /* return if stack is empty */
		    if (unlikely(curNode == BVH4i::invalidNode)) break;


		    /* intersect one ray against four triangles */
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
		      compactStack(stack_node_single,stack_dist_single,sindex,max_dist_xyz);

		  }	  
	      }
	    ray_tfar = select(valid0,ray16.tfar ,neg_inf);
	    continue;
	  }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	const unsigned int leaf_mask = BVH4I_LEAF_MASK;

	const mic_f time     = ray16.time;
	const mic_f one_time = (mic_f::one() - time);

        while (1)
        {
          /* test if this is a leaf node */
          if (unlikely(curNode.isLeaf(leaf_mask))) break;
          
          STAT3(normal.trav_nodes,1,popcnt(ray_tfar > curDist),16);
          const Node* __restrict__ const node = curNode.node(nodes);
          
          const BVH4mb::Node* __restrict__ const nodeMB = (BVH4mb::Node*)node;

          /* pop of next node */
          sptr_node--;
          sptr_dist--;
          curNode = *sptr_node; 
          curDist = *sptr_dist;
          
	  prefetch<PFHINT_L1>((char*)node + 0*64); 
	  prefetch<PFHINT_L1>((char*)node + 1*64); 
	  prefetch<PFHINT_L1>((char*)node + 2*64); 
	  prefetch<PFHINT_L1>((char*)node + 3*64); 

#pragma unroll(4)
          for (unsigned int i=0; i<4; i++)
          {
	    const NodeRef child = node->lower[i].child;

	    const mic_f lower_x =  one_time * nodeMB->lower[i].x + time * nodeMB->lower_t1[i].x;
	    const mic_f lower_y =  one_time * nodeMB->lower[i].y + time * nodeMB->lower_t1[i].y;
	    const mic_f lower_z =  one_time * nodeMB->lower[i].z + time * nodeMB->lower_t1[i].z;
	    const mic_f upper_x =  one_time * nodeMB->upper[i].x + time * nodeMB->upper_t1[i].x;
	    const mic_f upper_y =  one_time * nodeMB->upper[i].y + time * nodeMB->upper_t1[i].y;
	    const mic_f upper_z =  one_time * nodeMB->upper[i].z + time * nodeMB->upper_t1[i].z;

	    if (unlikely(i >=2 && child == BVH4i::invalidNode)) break;

            const mic_f lclipMinX = msub(lower_x,rdir16.x,org_rdir16.x);
            const mic_f lclipMinY = msub(lower_y,rdir16.y,org_rdir16.y);
            const mic_f lclipMinZ = msub(lower_z,rdir16.z,org_rdir16.z);
            const mic_f lclipMaxX = msub(upper_x,rdir16.x,org_rdir16.x);
            const mic_f lclipMaxY = msub(upper_y,rdir16.y,org_rdir16.y);
            const mic_f lclipMaxZ = msub(upper_z,rdir16.z,org_rdir16.z);
	    
            const mic_f lnearP = max(max(min(lclipMinX, lclipMaxX), min(lclipMinY, lclipMaxY)), min(lclipMinZ, lclipMaxZ));
            const mic_f lfarP  = min(min(max(lclipMinX, lclipMaxX), max(lclipMinY, lclipMaxY)), max(lclipMinZ, lclipMaxZ));
            const mic_m lhit   = max(lnearP,ray_tnear) <= min(lfarP,ray_tfar);   
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
              else 
		{
		  *(sptr_node-1) = child;
		  *(sptr_dist-1) = childDist; 
		}
              assert(sptr_node - stack_node < BVH4i::maxDepth);
            }	      
          }
#if SWITCH_ON_DOWN_TRAVERSAL == 1
	  const mic_m curUtil = ray_tfar > curDist;
	  if (unlikely(countbits(curUtil) <= BVH4i::hybridSIMDUtilSwitchThreshold))
	    {
	      *sptr_node++ = curNode;
	      *sptr_dist++ = curDist; 
	      goto pop;
	    }
#endif

        }
        
        /* return if stack is empty */
        if (unlikely(curNode == BVH4i::invalidNode)) break;
        

        /* intersect leaf */
        const mic_m m_valid_leaf = ray_tfar > curDist;
        STAT3(normal.trav_leaves,1,popcnt(m_valid_leaf),16);

	LeafIntersector::intersect16(curNode,
				     m_valid_leaf,
				     ray16.dir,
				     ray16.org,
				     ray16,
				     accel,
				     (Scene*)bvh->geometry);


        ray_tfar = select(m_valid_leaf,ray16.tfar,ray_tfar);
      }
    }

    template<typename LeafIntersector>    
    void BVH4mbIntersector16Hybrid<LeafIntersector>::occluded(mic_i* valid_i, BVH4mb* bvh, Ray16& ray16)
    {
      /* allocate stack */
      __aligned(64) mic_f   stack_dist[3*BVH4i::maxDepth+1];
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];
      __aligned(64) NodeRef stack_node_single[3*BVH4i::maxDepth+1];

      /* load ray */
      const mic_m m_valid     = *(mic_i*)valid_i != mic_i(0);
      mic_m m_terminated      = !m_valid;
      const mic3f rdir16      = rcp_safe(ray16.dir);
      const mic3f org_rdir16  = ray16.org * rdir16;
      mic_f ray_tnear         = select(m_valid,ray16.tnear,pos_inf);
      mic_f ray_tfar          = select(m_valid,ray16.tfar ,neg_inf);
      const mic_f inf = mic_f(pos_inf);

      
      /* push root node */
      stack_node[0] = BVH4i::invalidNode;
      stack_dist[0] = inf;
      stack_node[1] = bvh->root;
      stack_dist[1] = ray_tnear; 
      NodeRef* __restrict__ sptr_node = stack_node + 2;
      mic_f*   __restrict__ sptr_dist = stack_dist + 2;
      
      const Node      * __restrict__ nodes = (Node     *)bvh->nodePtr();
      const BVH4mb::Triangle01 * __restrict__ accel = (BVH4mb::Triangle01 *)bvh->triPtr();

      while (1) pop:
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

	/* switch to single ray mode */
        if (unlikely(countbits(m_stackDist) <= BVH4i::hybridSIMDUtilSwitchThreshold)) 
	  {
	    stack_node_single[0] = BVH4i::invalidNode;

	    /* traverse single ray */	  	  
	    long rayIndex = -1;
	    while((rayIndex = bitscan64(rayIndex,m_stackDist)) != BITSCAN_NO_BIT_SET_64) 
	      {	    
		stack_node_single[1] = curNode;
		size_t sindex = 2;

		const mic_f org_xyz      = loadAOS4to16f(rayIndex,ray16.org.x,ray16.org.y,ray16.org.z);
		const mic_f dir_xyz      = loadAOS4to16f(rayIndex,ray16.dir.x,ray16.dir.y,ray16.dir.z);
		const mic_f rdir_xyz     = loadAOS4to16f(rayIndex,rdir16.x,rdir16.y,rdir16.z);
		const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
		const mic_f min_dist_xyz = broadcast1to16f(&ray16.tnear[rayIndex]);
		const mic_f max_dist_xyz = broadcast1to16f(&ray16.tfar[rayIndex]);
		const mic_f time         = broadcast1to16f(&ray16.time[rayIndex]);

		const unsigned int leaf_mask = BVH4I_LEAF_MASK;

		while (1) 
		  {
		    NodeRef curNode = stack_node_single[sindex-1];
		    sindex--;
            
	  
		    traverse_single_occluded(curNode,
					     sindex,
					     rdir_xyz,
					     org_rdir_xyz,
					     min_dist_xyz,
					     max_dist_xyz,
					     time,
					     stack_node_single,
					     nodes,
					     leaf_mask);
	    

		    /* return if stack is empty */
		    if (unlikely(curNode == BVH4i::invalidNode)) break;

		    const mic_f zero = mic_f::zero();

		    /* intersect one ray against four triangles */

		    const bool hit = LeafIntersector::occluded(curNode,
							       rayIndex,
							       dir_xyz,
							       org_xyz,
							       min_dist_xyz,
							       max_dist_xyz,
							       ray16,
							       m_terminated,
							       accel,
							       (Scene*)bvh->geometry);
		    if (unlikely(hit)) break;
		  }	  

		if (unlikely(all(m_terminated))) 
		  {
		    store16i(m_valid,&ray16.geomID,mic_i::zero());
		    return;
		  }      

	      }
	    continue;
	  }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	const unsigned int leaf_mask = BVH4I_LEAF_MASK;

	const mic_f time     = ray16.time;
	const mic_f one_time = (mic_f::one() - time);


        while (1)
        {
          /* test if this is a leaf node */
          if (unlikely(curNode.isLeaf(leaf_mask))) break;
          
          STAT3(shadow.trav_nodes,1,popcnt(ray_tfar > curDist),16);
          const Node* __restrict__ const node = curNode.node(nodes);
          
          const BVH4mb::Node* __restrict__ const nodeMB = (BVH4mb::Node*)node;

	  prefetch<PFHINT_L1>((char*)node + 0*64);
	  prefetch<PFHINT_L1>((char*)node + 1*64);
	  prefetch<PFHINT_L1>((char*)node + 2*64);
	  prefetch<PFHINT_L1>((char*)node + 3*64);

          /* pop of next node */
          sptr_node--;
          sptr_dist--;
          curNode = *sptr_node; 
          curDist = *sptr_dist;
          
#pragma unroll(4)
          for (size_t i=0; i<4; i++)
	    {
	      const NodeRef child = node->lower[i].child;
            
	      const mic_f lower_x =  one_time * nodeMB->lower[i].x + time * nodeMB->lower_t1[i].x;
	      const mic_f lower_y =  one_time * nodeMB->lower[i].y + time * nodeMB->lower_t1[i].y;
	      const mic_f lower_z =  one_time * nodeMB->lower[i].z + time * nodeMB->lower_t1[i].z;
	      const mic_f upper_x =  one_time * nodeMB->upper[i].x + time * nodeMB->upper_t1[i].x;
	      const mic_f upper_y =  one_time * nodeMB->upper[i].y + time * nodeMB->upper_t1[i].y;
	      const mic_f upper_z =  one_time * nodeMB->upper[i].z + time * nodeMB->upper_t1[i].z;

	      if (unlikely(i >=2 && child == BVH4i::invalidNode)) break;

	      const mic_f lclipMinX = msub(lower_x,rdir16.x,org_rdir16.x);
	      const mic_f lclipMinY = msub(lower_y,rdir16.y,org_rdir16.y);
	      const mic_f lclipMinZ = msub(lower_z,rdir16.z,org_rdir16.z);
	      const mic_f lclipMaxX = msub(upper_x,rdir16.x,org_rdir16.x);
	      const mic_f lclipMaxY = msub(upper_y,rdir16.y,org_rdir16.y);
	      const mic_f lclipMaxZ = msub(upper_z,rdir16.z,org_rdir16.z);

	      const mic_f lnearP = max(max(min(lclipMinX, lclipMaxX), min(lclipMinY, lclipMaxY)), min(lclipMinZ, lclipMaxZ));
	      const mic_f lfarP  = min(min(max(lclipMinX, lclipMaxX), max(lclipMinY, lclipMaxY)), max(lclipMinZ, lclipMaxZ));
	      const mic_m lhit   = max(lnearP,ray_tnear) <= min(lfarP,ray_tfar);      
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
#if SWITCH_ON_DOWN_TRAVERSAL == 1
	  const mic_m curUtil = ray_tfar > curDist;
	  if (unlikely(countbits(curUtil) <= BVH4i::hybridSIMDUtilSwitchThreshold))
	    {
	      *sptr_node++ = curNode;
	      *sptr_dist++ = curDist; 
	      goto pop;
	    }
#endif
        }
        
        /* return if stack is empty */
        if (unlikely(curNode == BVH4i::invalidNode)) break;
        
        /* intersect leaf */
        mic_m m_valid_leaf = gt(m_active,ray_tfar,curDist);
        STAT3(shadow.trav_leaves,1,popcnt(m_valid_leaf),16);

	LeafIntersector::occluded16(curNode,
				    m_valid_leaf,
				    ray16.dir,
				    ray16.org,
				    ray16,
				    m_terminated,
				    accel,
				    (Scene*)bvh->geometry);

        ray_tfar = select(m_terminated,neg_inf,ray_tfar);
        if (unlikely(all(m_terminated))) break;
      }
      store16i(m_valid & m_terminated,&ray16.geomID,mic_i::zero());
    }
    
    DEFINE_INTERSECTOR16    (BVH4mbTriangle1Intersector16HybridMoeller, BVH4mbIntersector16Hybrid<Triangle1mbLeafIntersector>);
  }
}
