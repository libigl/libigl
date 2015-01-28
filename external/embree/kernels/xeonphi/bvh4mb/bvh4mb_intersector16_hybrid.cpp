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

#include "bvh4mb_intersector16_hybrid.h"
#include "geometry/triangle1.h"

namespace embree
{
  namespace isa
  {
    static unsigned int BVH4I_LEAF_MASK = BVH4i::leaf_mask; // needed due to compiler efficiency bug

    static __aligned(64) int zlc4[4] = {0xffffffff,0xffffffff,0xffffffff,0};

    void BVH4mbIntersector16Hybrid::intersect(mic_i* valid_i, BVH4mb* bvh, Ray16& ray16)
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
            
		    const mic_f one_time = (mic_f::one() - time);

		    while (1) 
		      {
			/* test if this is a leaf node */
			if (unlikely(curNode.isLeaf(leaf_mask))) break;
        
			const Node* __restrict__ const node = curNode.node(nodes);
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
			const mic_f tLowerXYZ = lower * rdir_xyz - org_rdir_xyz;
			const mic_f tUpperXYZ = upper * rdir_xyz - org_rdir_xyz;

			const mic_f tLower = mask_min(0x7777,min_dist_xyz,tLowerXYZ,tUpperXYZ);
			const mic_f tUpper = mask_max(0x7777,max_dist_xyz,tLowerXYZ,tUpperXYZ);

			sindex--;

			curNode = stack_node_single[sindex]; // early pop of next node

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
				stack_node_single[sindex] = node_second;
				((unsigned int*)stack_dist_single)[sindex] = dist_second;                      
				sindex++;
				assert(sindex < 3*BVH4i::maxDepth+1);
				continue;
			      }
			    else
			      {
				stack_node_single[sindex] = curNode;
				((unsigned int*)stack_dist_single)[sindex] = dist_first;
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
			const mic_i plower_node = load16i((int*)plower);
			const mic_m m_pos = andn(hitm,andn(closest_child,(mic_m)((unsigned int)closest_child - 1)));
			curNode = ((unsigned int*)plower)[closest_child_pos];

			compactustore16f(m_pos,&stack_dist_single[old_sindex],tNear); 
			compactustore16i(m_pos,&stack_node_single[old_sindex],plower_node);
		      }
	  
	    

		    /* return if stack is empty */
		    if (unlikely(curNode == BVH4i::invalidNode)) break;


		    /* intersect one ray against four triangles */
		    const mic_f zero = mic_f::zero();

		    const BVH4mb::Triangle01* tptr  = (BVH4mb::Triangle01*) curNode.leaf(accel);
	      
		    prefetch<PFHINT_L1>((mic_f*)tptr +  0); 
		    prefetch<PFHINT_L1>((mic_f*)tptr +  1); 
		    prefetch<PFHINT_L1>((mic_f*)tptr +  2); 
		    prefetch<PFHINT_L1>((mic_f*)tptr +  3); 
		    prefetch<PFHINT_L1>((mic_f*)tptr +  4); 
		    prefetch<PFHINT_L1>((mic_f*)tptr +  5); 
		    prefetch<PFHINT_L1>((mic_f*)tptr +  6); 
		    prefetch<PFHINT_L1>((mic_f*)tptr +  7); 

		    const mic_i and_mask = broadcast4to16i(zlc4);
	     

		    const mic_f v0_t0 = gather_4f_zlc(and_mask,
						      (float*)&tptr[0].t0.v0,
						      (float*)&tptr[1].t0.v0,
						      (float*)&tptr[2].t0.v0,
						      (float*)&tptr[3].t0.v0);
	      
		    const mic_f v1_t0 = gather_4f_zlc(and_mask,
						      (float*)&tptr[0].t0.v1,
						      (float*)&tptr[1].t0.v1,
						      (float*)&tptr[2].t0.v1,
						      (float*)&tptr[3].t0.v1);
	      
		    const mic_f v2_t0 = gather_4f_zlc(and_mask,
						      (float*)&tptr[0].t0.v2,
						      (float*)&tptr[1].t0.v2,
						      (float*)&tptr[2].t0.v2,
						      (float*)&tptr[3].t0.v2);

		    const mic_f v0_t1 = gather_4f_zlc(and_mask,
						      (float*)&tptr[0].t1.v0,
						      (float*)&tptr[1].t1.v0,
						      (float*)&tptr[2].t1.v0,
						      (float*)&tptr[3].t1.v0);
	      
		    const mic_f v1_t1 = gather_4f_zlc(and_mask,
						      (float*)&tptr[0].t1.v1,
						      (float*)&tptr[1].t1.v1,
						      (float*)&tptr[2].t1.v1,
						      (float*)&tptr[3].t1.v1);
	      
		    const mic_f v2_t1 = gather_4f_zlc(and_mask,
						      (float*)&tptr[0].t1.v2,
						      (float*)&tptr[1].t1.v2,
						      (float*)&tptr[2].t1.v2,
						      (float*)&tptr[3].t1.v2);

		    const mic_f v0 = v0_t0 * one_time + time * v0_t1;
		    const mic_f v1 = v1_t0 * one_time + time * v1_t1;
		    const mic_f v2 = v2_t0 * one_time + time * v2_t1;

		    const mic_f e1 = v1 - v0;
		    const mic_f e2 = v0 - v2;	     
		    const mic_f normal = lcross_zxy(e1,e2);
		    const mic_f org = v0 - org_xyz;
		    const mic_f odzxy = msubr231(org * swizzle(dir_xyz,_MM_SWIZ_REG_DACB), dir_xyz, swizzle(org,_MM_SWIZ_REG_DACB));
		    const mic_f den = ldot3_zxy(dir_xyz,normal);	      
		    const mic_f rcp_den = rcp(den);
		    const mic_f uu = ldot3_zxy(e2,odzxy); 
		    const mic_f vv = ldot3_zxy(e1,odzxy); 
		    const mic_f u = uu * rcp_den;
		    const mic_f v = vv * rcp_den;

#if defined(__BACKFACE_CULLING__)
		    const mic_m m_init = (mic_m)0x1111 & (den > zero);
#else
		    const mic_m m_init = 0x1111;
#endif

		    const mic_m valid_u = ge(m_init,u,zero);
		    const mic_m valid_v = ge(valid_u,v,zero);
		    const mic_m m_aperture = le(valid_v,u+v,mic_f::one()); 

		    const mic_f nom = ldot3_zxy(org,normal);
		    if (unlikely(none(m_aperture))) continue;
		    const mic_f t = rcp_den*nom;

		    mic_m m_final  = lt(lt(m_aperture,min_dist_xyz,t),t,max_dist_xyz);

		    max_dist_xyz  = select(m_final,t,max_dist_xyz);		    

#if defined(__USE_RAY_MASK__)
		    const mic_i rayMask(ray16.mask[rayIndex]);
		    const mic_i triMask = swDDDD(gather16i_4i_align(&tptr[0].t0.v2,&tptr[1].t0.v2,&tptr[2].t0.v2,&tptr[3].t0.v2));
		    const mic_m m_ray_mask = (rayMask & triMask) != mic_i::zero();
		    m_final &= m_ray_mask;	      
#endif

		    /* did the ray hot one of the four triangles? */
		    if (unlikely(any(m_final)))
		      {
			prefetch<PFHINT_L1EX>(&ray16.tfar);  
			prefetch<PFHINT_L1EX>(&ray16.u);
			prefetch<PFHINT_L1EX>(&ray16.v);
			prefetch<PFHINT_L1EX>(&ray16.Ng.x); 
			prefetch<PFHINT_L1EX>(&ray16.Ng.y); 
			prefetch<PFHINT_L1EX>(&ray16.Ng.z); 
			prefetch<PFHINT_L1EX>(&ray16.geomID);
			prefetch<PFHINT_L1EX>(&ray16.primID);

			const mic_f min_dist = vreduce_min(max_dist_xyz);
			const mic_m m_dist = eq(min_dist,max_dist_xyz);

			const size_t vecIndex = bitscan(toInt(m_dist));
			const size_t triIndex = vecIndex >> 2;

			const BVH4mb::Triangle01  *__restrict__ tri_ptr = tptr + triIndex;

			const mic_m m_tri = m_dist^(m_dist & (mic_m)((unsigned int)m_dist - 1));

			const mic_f gnormalz = swAAAA(normal);
			const mic_f gnormalx = swBBBB(normal);
			const mic_f gnormaly = swCCCC(normal);
		  
			max_dist_xyz = min_dist;

			compactustore16f_low(m_tri,&ray16.tfar[rayIndex],min_dist);
			compactustore16f_low(m_tri,&ray16.u[rayIndex],u); 
			compactustore16f_low(m_tri,&ray16.v[rayIndex],v); 
			compactustore16f_low(m_tri,&ray16.Ng.x[rayIndex],gnormalx); 
			compactustore16f_low(m_tri,&ray16.Ng.y[rayIndex],gnormaly); 
			compactustore16f_low(m_tri,&ray16.Ng.z[rayIndex],gnormalz); 

			ray16.geomID[rayIndex] = tri_ptr->t0.geomID();
			ray16.primID[rayIndex] = tri_ptr->t0.primID();

			/* compact the stack if size of stack >= 2 */
			if (likely(sindex >= 2))
			  {
			    if (likely(sindex < 16))
			      {
				const unsigned int m_num_stack = mic_m::shift1[sindex] - 1;
				const mic_m m_num_stack_low  = toMask(m_num_stack);
				const mic_f snear_low  = load16f(stack_dist_single + 0);
				const mic_i snode_low  = load16i((int*)stack_node_single + 0);
				const mic_m m_stack_compact_low  = le(m_num_stack_low,snear_low,max_dist_xyz) | (mic_m)1;
				compactustore16f_low(m_stack_compact_low,stack_dist_single + 0,snear_low);
				compactustore16i_low(m_stack_compact_low,(int*)stack_node_single + 0,snode_low);
				sindex = countbits(m_stack_compact_low);
				assert(sindex < 16);
			      }
			    else if (likely(sindex < 32))
			      {
				const mic_m m_num_stack_high = toMask(mic_m::shift1[sindex-16] - 1); 
				const mic_f snear_low  = load16f(stack_dist_single + 0);
				const mic_f snear_high = load16f(stack_dist_single + 16);
				const mic_i snode_low  = load16i((int*)stack_node_single + 0);
				const mic_i snode_high = load16i((int*)stack_node_single + 16);
				const mic_m m_stack_compact_low  = le(snear_low,max_dist_xyz) | (mic_m)1;
				const mic_m m_stack_compact_high = le(m_num_stack_high,snear_high,max_dist_xyz);
				compactustore16f(m_stack_compact_low,      stack_dist_single + 0,snear_low);
				compactustore16i(m_stack_compact_low,(int*)stack_node_single + 0,snode_low);
				compactustore16f(m_stack_compact_high,      stack_dist_single + countbits(m_stack_compact_low),snear_high);
				compactustore16i(m_stack_compact_high,(int*)stack_node_single + countbits(m_stack_compact_low),snode_high);
				assert ((unsigned int)m_num_stack_high == ((mic_m::shift1[sindex] - 1) >> 16));
				sindex = countbits(m_stack_compact_low) + countbits(m_stack_compact_high);
				assert(sindex < 32);
			      }
			    else
			      {
				const mic_m m_num_stack_32 = toMask(mic_m::shift1[sindex-32] - 1); 

				const mic_f snear_0  = load16f(stack_dist_single + 0);
				const mic_f snear_16 = load16f(stack_dist_single + 16);
				const mic_f snear_32 = load16f(stack_dist_single + 32);
				const mic_i snode_0  = load16i((int*)stack_node_single + 0);
				const mic_i snode_16 = load16i((int*)stack_node_single + 16);
				const mic_i snode_32 = load16i((int*)stack_node_single + 32);
				const mic_m m_stack_compact_0  = le(               snear_0 ,max_dist_xyz) | (mic_m)1;
				const mic_m m_stack_compact_16 = le(               snear_16,max_dist_xyz);
				const mic_m m_stack_compact_32 = le(m_num_stack_32,snear_32,max_dist_xyz);

				sindex = 0;
				compactustore16f(m_stack_compact_0,      stack_dist_single + sindex,snear_0);
				compactustore16i(m_stack_compact_0,(int*)stack_node_single + sindex,snode_0);
				sindex += countbits(m_stack_compact_0);
				compactustore16f(m_stack_compact_16,      stack_dist_single + sindex,snear_16);
				compactustore16i(m_stack_compact_16,(int*)stack_node_single + sindex,snode_16);
				sindex += countbits(m_stack_compact_16);
				compactustore16f(m_stack_compact_32,      stack_dist_single + sindex,snear_32);
				compactustore16i(m_stack_compact_32,(int*)stack_node_single + sindex,snode_32);
				sindex += countbits(m_stack_compact_32);

				assert(sindex < 48);		  
			      }
			  }
		      }
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

            //if (unlikely(child == BVH4i::emptyNode)) break;

	    const mic_f lower_x =  one_time * nodeMB->lower[i].x + time * nodeMB->lower_t1[i].x;
	    const mic_f lower_y =  one_time * nodeMB->lower[i].y + time * nodeMB->lower_t1[i].y;
	    const mic_f lower_z =  one_time * nodeMB->lower[i].z + time * nodeMB->lower_t1[i].z;
	    const mic_f upper_x =  one_time * nodeMB->upper[i].x + time * nodeMB->upper_t1[i].x;
	    const mic_f upper_y =  one_time * nodeMB->upper[i].y + time * nodeMB->upper_t1[i].y;
	    const mic_f upper_z =  one_time * nodeMB->upper[i].z + time * nodeMB->upper_t1[i].z;


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
        }
        
        /* return if stack is empty */
        if (unlikely(curNode == BVH4i::invalidNode)) break;
        
	const mic3f org = ray16.org;

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
	    const mic_f den = dot(Ng,ray16.dir);

	    const mic_f rcp_den = rcp(den);

	    mic_m valid = valid_leaf;

#if defined(__BACKFACE_CULLING__)
	    
	    valid &= den > zero;
#endif

	    /* perform edge tests */
	    const mic3f R = -cross(C,ray16.dir);
	    const mic_f u = dot(R,e2)*rcp_den;
	    const mic_f v = dot(R,e1)*rcp_den;
	    valid = ge(valid,u,zero);
	    valid = ge(valid,v,zero);
	    valid = le(valid,u+v,one);

	    prefetch<PFHINT_L1EX>(&ray16.u);      
	    prefetch<PFHINT_L1EX>(&ray16.v);      
	    prefetch<PFHINT_L1EX>(&ray16.tfar);      


	    if (unlikely(none(valid))) continue;

	    const mic_f dot_C_Ng = dot(C,Ng);
	    const mic_f t = dot_C_Ng * rcp_den;
      
	    /* perform depth test */
	    valid = ge(valid, t,ray16.tnear);
	    valid = ge(valid,ray16.tfar,t);

	    const mic_i geomID = tri_t0.geomID();
	    const mic_i primID = tri_t0.primID();
	    prefetch<PFHINT_L1EX>(&ray16.geomID);      
	    prefetch<PFHINT_L1EX>(&ray16.primID);      
	    prefetch<PFHINT_L1EX>(&ray16.Ng.x);      
	    prefetch<PFHINT_L1EX>(&ray16.Ng.y);      
	    prefetch<PFHINT_L1EX>(&ray16.Ng.z);      

	    /* ray masking test */
#if defined(__USE_RAY_MASK__)
	    valid &= (tri_t0.mask() & ray16.mask) != 0;
#endif
	    if (unlikely(none(valid))) continue;

	    /* update hit information */
	    store16f(valid,(float*)&ray16.u,u);
	    store16f(valid,(float*)&ray16.v,v);
	    store16f(valid,(float*)&ray16.tfar,t);
	    store16i(valid,(float*)&ray16.geomID,geomID);
	    store16i(valid,(float*)&ray16.primID,primID);
	    store16f(valid,(float*)&ray16.Ng.x,Ng.x);
	    store16f(valid,(float*)&ray16.Ng.y,Ng.y);
	    store16f(valid,(float*)&ray16.Ng.z,Ng.z);
	  }

        ray_tfar = select(valid_leaf,ray16.tfar,ray_tfar);
      }
    }
    
    void BVH4mbIntersector16Hybrid::occluded(mic_i* valid_i, BVH4mb* bvh, Ray16& ray16)
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
            
		    const mic_f one_time = (mic_f::one() - time);

		    while (1) 
		      {
			/* test if this is a leaf node */
			if (unlikely(curNode.isLeaf(leaf_mask))) break;
        
			const Node* __restrict__ const node = curNode.node(nodes);
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
			const mic_f tLowerXYZ = lower * rdir_xyz - org_rdir_xyz;
			const mic_f tUpperXYZ = upper * rdir_xyz - org_rdir_xyz;

			const mic_f tLower = mask_min(0x7777,min_dist_xyz,tLowerXYZ,tUpperXYZ);
			const mic_f tUpper = mask_max(0x7777,max_dist_xyz,tLowerXYZ,tUpperXYZ);

			sindex--;
			curNode = stack_node_single[sindex]; // early pop of next node

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
				stack_node_single[sindex] = node_second;
				sindex++;
				assert(sindex < 3*BVH4i::maxDepth+1);
				continue;
			      }
			    else
			      {
				stack_node_single[sindex] = curNode;
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
			compactustore16i(m_pos,&stack_node_single[old_sindex],plower_node);
		      }
	  
	    

		    /* return if stack is empty */
		    if (unlikely(curNode == BVH4i::invalidNode)) break;

		    const mic_f zero = mic_f::zero();

		    /* intersect one ray against four triangles */

		    const BVH4mb::Triangle01* tptr  = (BVH4mb::Triangle01*) curNode.leaf(accel);

		    prefetch<PFHINT_L1>((mic_f*)tptr +  0); 
		    prefetch<PFHINT_L1>((mic_f*)tptr +  1); 
		    prefetch<PFHINT_L1>((mic_f*)tptr +  2); 
		    prefetch<PFHINT_L1>((mic_f*)tptr +  3); 
		    prefetch<PFHINT_L1>((mic_f*)tptr +  4); 
		    prefetch<PFHINT_L1>((mic_f*)tptr +  5); 
		    prefetch<PFHINT_L1>((mic_f*)tptr +  6); 
		    prefetch<PFHINT_L1>((mic_f*)tptr +  7); 

		    const mic_i and_mask = broadcast4to16i(zlc4);
	      
		    const mic_f v0_t0 = gather_4f_zlc(and_mask,
						      (float*)&tptr[0].t0.v0,
						      (float*)&tptr[1].t0.v0,
						      (float*)&tptr[2].t0.v0,
						      (float*)&tptr[3].t0.v0);
	      
		    const mic_f v1_t0 = gather_4f_zlc(and_mask,
						      (float*)&tptr[0].t0.v1,
						      (float*)&tptr[1].t0.v1,
						      (float*)&tptr[2].t0.v1,
						      (float*)&tptr[3].t0.v1);
	      
		    const mic_f v2_t0 = gather_4f_zlc(and_mask,
						      (float*)&tptr[0].t0.v2,
						      (float*)&tptr[1].t0.v2,
						      (float*)&tptr[2].t0.v2,
						      (float*)&tptr[3].t0.v2);

		    const mic_f v0_t1 = gather_4f_zlc(and_mask,
						      (float*)&tptr[0].t1.v0,
						      (float*)&tptr[1].t1.v0,
						      (float*)&tptr[2].t1.v0,
						      (float*)&tptr[3].t1.v0);
	      
		    const mic_f v1_t1 = gather_4f_zlc(and_mask,
						      (float*)&tptr[0].t1.v1,
						      (float*)&tptr[1].t1.v1,
						      (float*)&tptr[2].t1.v1,
						      (float*)&tptr[3].t1.v1);
	      
		    const mic_f v2_t1 = gather_4f_zlc(and_mask,
						      (float*)&tptr[0].t1.v2,
						      (float*)&tptr[1].t1.v2,
						      (float*)&tptr[2].t1.v2,
						      (float*)&tptr[3].t1.v2);

		    const mic_f v0 = v0_t0 * one_time + time * v0_t1;
		    const mic_f v1 = v1_t0 * one_time + time * v1_t1;
		    const mic_f v2 = v2_t0 * one_time + time * v2_t1;

		    const mic_f e1 = v1 - v0;
		    const mic_f e2 = v0 - v2;	     
		    const mic_f normal = lcross_zxy(e1,e2);
		    const mic_f org = v0 - org_xyz;
		    const mic_f odzxy = msubr231(org * swizzle(dir_xyz,_MM_SWIZ_REG_DACB), dir_xyz, swizzle(org,_MM_SWIZ_REG_DACB));
		    const mic_f den = ldot3_zxy(dir_xyz,normal);	      
		    const mic_f rcp_den = rcp(den);
		    const mic_f uu = ldot3_zxy(e2,odzxy); 
		    const mic_f vv = ldot3_zxy(e1,odzxy); 
		    const mic_f u = uu * rcp_den;
		    const mic_f v = vv * rcp_den;

#if defined(__BACKFACE_CULLING__)
		    const mic_m m_init = (mic_m)0x1111 & (den > zero);
#else
		    const mic_m m_init = 0x1111;
#endif

		    const mic_m valid_u = ge(m_init,u,zero);
		    const mic_m valid_v = ge(valid_u,v,zero);
		    const mic_m m_aperture = le(valid_v,u+v,mic_f::one()); 

		    const mic_f nom = ldot3_zxy(org,normal);
		    const mic_f t = rcp_den*nom;

		    if (unlikely(none(m_aperture))) continue;

		    mic_m m_final  = lt(lt(m_aperture,min_dist_xyz,t),t,max_dist_xyz);

#if defined(__USE_RAY_MASK__)
		    const mic_i rayMask(ray16.mask[rayIndex]);
		    const mic_i triMask = swDDDD(gather16i_4i_align(&tptr[0].t0.v2,&tptr[1].t0.v2,&tptr[2].t0.v2,&tptr[3].t0.v2));
		    const mic_m m_ray_mask = (rayMask & triMask) != mic_i::zero();
		    m_final &= m_ray_mask;	      
#endif

		    /* did the ray hot one of the four triangles? */
		    if (unlikely(any(m_final)))
		      {
			m_terminated |= toMask(mic_m::shift1[rayIndex]);
			break;
		      }
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

	mic_m valid0 = valid_leaf;

	const mic3f org = ray16.org;
	const mic3f dir = ray16.dir;

	const mic_f zero = mic_f::zero();
	const mic_f one  = mic_f::one();

     
	for (size_t i=0; i<items; i++) 
	  {
	    STAT3(shadow.trav_prims,1,popcnt(valid0),16);

	    const Triangle1& tri_t0 = tris[i].t0;
	    const Triangle1& tri_t1 = tris[i].t1;

	    prefetch<PFHINT_L1>(&tris[i+1].t0); 
	    prefetch<PFHINT_L1>(&tris[i+1].t1); 

	    mic_m valid = valid0;
        
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

	    const mic_f den = dot(dir,Ng);

#if defined(__BACKFACE_CULLING__)
	    valid &= den > zero;
#endif
	    const mic_f rcp_den = rcp(den);
	    const mic3f R = cross(dir,C);
	    const mic_f u = dot(R,e1)*rcp_den;
	    const mic_f v = dot(R,e2)*rcp_den;
	    valid = ge(valid,u,zero);
	    valid = ge(valid,v,zero);
	    valid = le(valid,u+v,one); 
	    const mic_f t = dot(C,Ng) * rcp_den;
	    evictL1(tris);

	    if (unlikely(none(valid))) continue;
      
	    /* perform depth test */
	    valid = ge(valid, t,ray16.tnear);
	    valid = ge(valid,ray16.tfar,t);

	    /* ray masking test */
#if defined(__USE_RAY_MASK__)
	    valid &= (tri_t0.mask() & ray16.mask) != 0;
#endif
	    if (unlikely(none(valid))) continue;

	    /* update occlusion */
	    valid0 &= !valid;
	    if (unlikely(none(valid0))) break;
	  }
	m_terminated |= valid_leaf & (!valid0);	

        ray_tfar = select(m_terminated,neg_inf,ray_tfar);
        if (unlikely(all(m_terminated))) break;
      }
      store16i(m_valid & m_terminated,&ray16.geomID,mic_i::zero());
    }
    
    DEFINE_INTERSECTOR16    (BVH4mbTriangle1Intersector16HybridMoeller, BVH4mbIntersector16Hybrid);
  }
}
