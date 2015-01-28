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

#include "bvh4mb_intersector16_chunk.h"
#include "geometry/triangle1.h"

namespace embree
{
  namespace isa
  {
    static unsigned int BVH4I_LEAF_MASK = BVH4i::leaf_mask; // needed due to compiler efficiency bug

    void BVH4mbIntersector16Chunk::intersect(mic_i* valid_i, BVH4mb* bvh, Ray16& ray)
    {
      /* near and node stack */
      __aligned(64) mic_f   stack_dist[3*BVH4i::maxDepth+1];
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

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
      
      const Node               * __restrict__ nodes = (Node     *)bvh->nodePtr();
      const BVH4mb::Triangle01 * __restrict__ accel = (BVH4mb::Triangle01 *)bvh->triPtr();

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
        if (unlikely(none(m_stackDist))) {continue;}
	        
	const unsigned int leaf_mask = BVH4I_LEAF_MASK; 

	const mic_f time     = ray.time;
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

	  prefetch<PFHINT_L1>((mic_f*)node + 0);           
	  prefetch<PFHINT_L1>((mic_f*)node + 1); 
	  prefetch<PFHINT_L1>((mic_f*)node + 2); 
	  prefetch<PFHINT_L1>((mic_f*)node + 3); 

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
              assert(sptr_node - stack_node < BVH4i::maxDepth);
            }	      
          }
        }
        
        /* return if stack is empty */
        if (unlikely(curNode == BVH4i::invalidNode)) break;
        
        /* intersect leaf */
        const mic_m valid_leaf = ray_tfar > curDist;
        STAT3(normal.trav_leaves,1,popcnt(valid_leaf),16);
 
	unsigned int items; 
	const BVH4mb::Triangle01* tris  = (BVH4mb::Triangle01*) curNode.leaf(accel,items);

	const mic_f zero = mic_f::zero();
	const mic_f one  = mic_f::one();

	prefetch<PFHINT_L1>((mic_f*)tris +  0); 
	prefetch<PFHINT_L2>((mic_f*)tris +  1); 
	prefetch<PFHINT_L2>((mic_f*)tris +  2); 
	prefetch<PFHINT_L2>((mic_f*)tris +  3); 
	prefetch<PFHINT_L2>((mic_f*)tris +  4); 
	prefetch<PFHINT_L2>((mic_f*)tris +  5); 
	prefetch<PFHINT_L2>((mic_f*)tris +  6); 
	prefetch<PFHINT_L2>((mic_f*)tris +  7); 

        const mic3f org = ray.org;
        const mic3f dir = ray.dir;

	for (size_t i=0; i<items; i++) 
	  {
	    const Triangle1& tri_t0 = tris[i].t0;
	    const Triangle1& tri_t1 = tris[i].t1;

	    prefetch<PFHINT_L1>(&tris[i+1].t0); 
	    prefetch<PFHINT_L1>(&tris[i+1].t1); 

	    STAT3(normal.trav_prims,1,popcnt(valid_i),16);
        
	    /* load vertices and calculate edges */
	    const mic3f v0_t0( broadcast1to16f(&tri_t0.v0.x), broadcast1to16f(&tri_t0.v0.y), broadcast1to16f(&tri_t0.v0.z) );
	    const mic3f v0_t1( broadcast1to16f(&tri_t1.v0.x), broadcast1to16f(&tri_t1.v0.y), broadcast1to16f(&tri_t1.v0.z) );
	    const mic3f v0 = v0_t0 * one_time + time * v0_t1;
	    const mic3f v1_t0( broadcast1to16f(&tri_t0.v1.x), broadcast1to16f(&tri_t0.v1.y), broadcast1to16f(&tri_t0.v1.z) );
	    const mic3f v1_t1( broadcast1to16f(&tri_t1.v1.x), broadcast1to16f(&tri_t1.v1.y), broadcast1to16f(&tri_t1.v1.z) );
	    const mic3f v1 = v1_t0 * one_time + time * v1_t1;
	    const mic3f v2_t0( broadcast1to16f(&tri_t0.v2.x), broadcast1to16f(&tri_t0.v2.y), broadcast1to16f(&tri_t0.v2.z) );
	    const mic3f v2_t1( broadcast1to16f(&tri_t1.v2.x), broadcast1to16f(&tri_t1.v2.y), broadcast1to16f(&tri_t1.v2.z) );
	    const mic3f v2 = v2_t0 * one_time + time * v2_t1;

	    const mic3f e1 = v0-v1;
	    const mic3f e2 = v2-v0;

	    const mic3f Ng = cross(e1,e2);

	    /* calculate denominator */
	    const mic3f C =  v0 - org;
	    
	    const mic_f den = dot(Ng,dir);

	    mic_m valid = valid_leaf;

#if defined(__BACKFACE_CULLING__)
	    
	    valid &= den > zero;
#endif

	    /* perform edge tests */
	    const mic_f rcp_den = rcp(den);
	    const mic3f R = cross(dir,C);
	    const mic_f u = dot(R,e2)*rcp_den;
	    const mic_f v = dot(R,e1)*rcp_den;
	    valid = ge(valid,u,zero);
	    valid = ge(valid,v,zero);
	    valid = le(valid,u+v,one);
	    prefetch<PFHINT_L1EX>(&ray.u);      
	    prefetch<PFHINT_L1EX>(&ray.v);      
	    prefetch<PFHINT_L1EX>(&ray.tfar);      
	    const mic_f t = dot(C,Ng) * rcp_den;

	    if (unlikely(none(valid))) continue;
      
	    /* perform depth test */
	    valid = ge(valid, t,ray.tnear);
	    valid = ge(valid,ray.tfar,t);

	    const mic_i geomID = tri_t0.geomID();
	    const mic_i primID = tri_t0.primID();
	    prefetch<PFHINT_L1EX>(&ray.geomID);      
	    prefetch<PFHINT_L1EX>(&ray.primID);      
	    prefetch<PFHINT_L1EX>(&ray.Ng.x);      
	    prefetch<PFHINT_L1EX>(&ray.Ng.y);      
	    prefetch<PFHINT_L1EX>(&ray.Ng.z);      

	    /* ray masking test */
#if defined(__USE_RAY_MASK__)
	    valid &= (mic_i(tri_t0.mask()) & ray.mask) != 0;
#endif
	    if (unlikely(none(valid))) continue;
        
	    /* update hit information */
	    store16f(valid,(float*)&ray.u,u);
	    store16f(valid,(float*)&ray.v,v);
	    store16f(valid,(float*)&ray.tfar,t);
	    store16i(valid,(float*)&ray.geomID,geomID);
	    store16i(valid,(float*)&ray.primID,primID);
	    store16f(valid,(float*)&ray.Ng.x,Ng.x);
	    store16f(valid,(float*)&ray.Ng.y,Ng.y);
	    store16f(valid,(float*)&ray.Ng.z,Ng.z);
	  }
        ray_tfar = select(valid_leaf,ray.tfar,ray_tfar);
      }
    }

    void BVH4mbIntersector16Chunk::occluded(mic_i* valid_i, BVH4mb* bvh, Ray16& ray)
    {
      /* allocate stack */
      __aligned(64) mic_f    stack_dist[3*BVH4i::maxDepth+1];
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

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
      
      const Node               * __restrict__ nodes = (Node     *)bvh->nodePtr();
      const BVH4mb::Triangle01 * __restrict__ accel = (BVH4mb::Triangle01 *)bvh->triPtr();

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

        if (unlikely(none(m_stackDist))) { continue; }
	
	const unsigned int leaf_mask = BVH4I_LEAF_MASK; 

	const mic_f time     = ray.time;
	const mic_f one_time = (mic_f::one() - time);

        while (1)
        {
          /* test if this is a leaf node */
          if (unlikely(curNode.isLeaf(leaf_mask))) break;
          
          STAT3(shadow.trav_nodes,1,popcnt(ray_tfar > curDist),16);
          const Node* __restrict__ const node = curNode.node(nodes);
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
	    const NodeRef child = node->lower[i].child;

	    const mic_f lower_x =  one_time * nodeMB->lower[i].x + time * nodeMB->lower_t1[i].x;
	    const mic_f lower_y =  one_time * nodeMB->lower[i].y + time * nodeMB->lower_t1[i].y;
	    const mic_f lower_z =  one_time * nodeMB->lower[i].z + time * nodeMB->lower_t1[i].z;
	    const mic_f upper_x =  one_time * nodeMB->upper[i].x + time * nodeMB->upper_t1[i].x;
	    const mic_f upper_y =  one_time * nodeMB->upper[i].y + time * nodeMB->upper_t1[i].y;
	    const mic_f upper_z =  one_time * nodeMB->upper[i].z + time * nodeMB->upper_t1[i].z;


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
              assert(sptr_node - stack_node < BVH4i::maxDepth);
            }	      
          }
        }
        
        /* return if stack is empty */
        if (unlikely(curNode == BVH4i::invalidNode)) break;
        
        /* intersect leaf */
        mic_m valid_leaf = gt(m_active,ray_tfar,curDist);
        STAT3(shadow.trav_leaves,1,popcnt(valid_leaf),16);

	unsigned int items; 
	const BVH4mb::Triangle01* tris  = (BVH4mb::Triangle01*) curNode.leaf(accel,items);

	prefetch<PFHINT_L1>((mic_f*)tris +  0); 
	prefetch<PFHINT_L2>((mic_f*)tris +  1); 
	prefetch<PFHINT_L2>((mic_f*)tris +  2); 
	prefetch<PFHINT_L2>((mic_f*)tris +  3); 
	prefetch<PFHINT_L2>((mic_f*)tris +  4); 
	prefetch<PFHINT_L2>((mic_f*)tris +  5); 
	prefetch<PFHINT_L2>((mic_f*)tris +  6); 
	prefetch<PFHINT_L2>((mic_f*)tris +  7); 

        const mic3f org = ray.org;
        const mic3f dir = ray.dir;
	const mic_f zero = mic_f::zero();

	for (size_t i=0; i<items; i++) 
	  {
	    const Triangle1& tri_t0 = tris[i].t0;
	    const Triangle1& tri_t1 = tris[i].t1;

	    prefetch<PFHINT_L1>(&tris[i+1].t0); 
	    prefetch<PFHINT_L1>(&tris[i+1].t1); 

	    STAT3(normal.trav_prims,1,popcnt(valid_i),16);
        
	    /* load vertices and calculate edges */
	    const mic3f v0_t0( broadcast1to16f(&tri_t0.v0.x), broadcast1to16f(&tri_t0.v0.y), broadcast1to16f(&tri_t0.v0.z) );
	    const mic3f v0_t1( broadcast1to16f(&tri_t1.v0.x), broadcast1to16f(&tri_t1.v0.y), broadcast1to16f(&tri_t1.v0.z) );
	    const mic3f v0 = v0_t0 * one_time + time * v0_t1;
	    const mic3f v1_t0( broadcast1to16f(&tri_t0.v1.x), broadcast1to16f(&tri_t0.v1.y), broadcast1to16f(&tri_t0.v1.z) );
	    const mic3f v1_t1( broadcast1to16f(&tri_t1.v1.x), broadcast1to16f(&tri_t1.v1.y), broadcast1to16f(&tri_t1.v1.z) );
	    const mic3f v1 = v1_t0 * one_time + time * v1_t1;
	    const mic3f v2_t0( broadcast1to16f(&tri_t0.v2.x), broadcast1to16f(&tri_t0.v2.y), broadcast1to16f(&tri_t0.v2.z) );
	    const mic3f v2_t1( broadcast1to16f(&tri_t1.v2.x), broadcast1to16f(&tri_t1.v2.y), broadcast1to16f(&tri_t1.v2.z) );
	    const mic3f v2 = v2_t0 * one_time + time * v2_t1;

	    const mic3f e1 = v0-v1;
	    const mic3f e2 = v2-v0;

	    const mic3f Ng = cross(e1,e2);

	    /* calculate denominator */
	    const mic3f C =  v0 - org;
	    
	    const mic_f den = dot(Ng,dir);

	    mic_m valid = valid_leaf;

#if defined(__BACKFACE_CULLING__)
	    
	    valid &= den > zero;
#endif

	    /* perform edge tests */
	    const mic_f rcp_den = rcp(den);
	    const mic3f R = cross(dir,C);
	    const mic_f u = dot(R,e2)*rcp_den;
	    const mic_f v = dot(R,e1)*rcp_den;
	    valid = ge(valid,u,zero);
	    valid = ge(valid,v,zero);
	    valid = le(valid,u+v,one);
	    const mic_f t = dot(C,Ng) * rcp_den;

	    if (unlikely(none(valid))) continue;
      
	    /* perform depth test */
	    valid = ge(valid, t,ray.tnear);
	    valid = ge(valid,ray.tfar,t);

	    /* ray masking test */
#if defined(__USE_RAY_MASK__)
	    valid &= (mic_i(tri_t0.mask()) & ray.mask) != 0;
#endif
	    if (unlikely(none(valid))) continue;
	    
	    /* update occlusion */
	    m_terminated |= valid;
	    valid_leaf &= ~valid;
	    if (unlikely(none(valid_leaf))) break;
	  }
        if (unlikely(all(m_terminated))) break;
        ray_tfar = select(m_terminated,neg_inf,ray_tfar);
      }
      store16i(valid & m_terminated,&ray.geomID,0);
    }

    DEFINE_INTERSECTOR16    (BVH4mbTriangle1Intersector16ChunkMoeller, BVH4mbIntersector16Chunk);
  }
}
