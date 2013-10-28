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

#include "bvh4aos_triangle1_intersector16_hybrid_ref.h"

#define STAT_COLLECTOR(x)

#define SIMD_UTIL_SWITCH_THRESHOLD 8
#define QUAD_BVH_MAX_STACK_DEPTH 64

#define BVH_INDEX_SHIFT 6
#define BVH_ITEMS_MASK   (((unsigned int)1 << BVH_INDEX_SHIFT)-1)
#define BVH_LEAF_MASK    ((unsigned int)1 << 31)
#define BVH_OFFSET_MASK  (~(BVH_ITEMS_MASK | BVH_LEAF_MASK))

#define QBVH_INDEX_SHIFT 7
#define QBVH_LEAF_BIT_SHIFT 5
#define QBVH_ITEMS_MASK   (((unsigned int)1 << QBVH_LEAF_BIT_SHIFT)-1)
#define QBVH_LEAF_MASK     ((unsigned int)1 << QBVH_LEAF_BIT_SHIFT)
#define QBVH_OFFSET_MASK  (~(QBVH_ITEMS_MASK | QBVH_LEAF_MASK))
#define QBVH_TERMINAL_TOKEN QBVH_LEAF_MASK

//#define USE_PRECOMPUTED_GNORMAL

namespace embree
{
  __forceinline unsigned int qbvhItemOffset(const unsigned int children) {
    return children & BVH_OFFSET_MASK; // 6 bits instead of 7
  }
  
  __forceinline unsigned int qbvhItemOffsetToID(const unsigned int children) {
    return children >> BVH_INDEX_SHIFT; // 6 bits instead of 7
  }
  
  __forceinline unsigned int qbvhItems(const unsigned int children) {
    return children & QBVH_ITEMS_MASK; // 6 bits instead of 7
  }
  
  __forceinline unsigned int qbvhChildID(const unsigned int node) {
    return (node & QBVH_OFFSET_MASK) >> QBVH_INDEX_SHIFT;
  }

  __forceinline BVH4AOS::Node *qbvhChildPtr(const BVH4AOS::Node * __restrict__ const ptr, const unsigned int node) {
    const unsigned int offset = node & QBVH_OFFSET_MASK;
    return (BVH4AOS::Node*)((char*)ptr + offset);
  }

  __forceinline BVH4AOS::Node *qbvhChildPtrNoMask(const BVH4AOS::Node * __restrict__ const ptr, const unsigned int node) {
    return (BVH4AOS::Node*)((char*)ptr + (unsigned long)node);
  }

  __forceinline unsigned int qbvhLeaf(const unsigned int node) {
    return (node & QBVH_LEAF_MASK);
  }

  __forceinline unsigned int qbvhLeaf(const unsigned int node, const unsigned int mask) {
    return (node & mask);
  }

  __forceinline unsigned int qbvhChildren(const unsigned int node) {
    return (node & QBVH_ITEMS_MASK);
  }
  
  __forceinline unsigned int qbvhCreateNode(const unsigned int nodeID, const unsigned int children) {
    return (nodeID << QBVH_INDEX_SHIFT) | children;
  };

  static unsigned int TMP_QBVH_LEAF_MASK = /*TMP_*/QBVH_LEAF_MASK; // needed due to compiler efficiency bug

  void BVH4AOSTriangle1Intersector16HybridRef::intersect(const BVH4AOSTriangle1Intersector16HybridRef* This, Ray16& ray, const __mmask m_active_i) 
  {
    MIC_ALIGN float        pack_near[16];
    MIC_ALIGN unsigned int stack_node_chunk[QUAD_BVH_MAX_STACK_DEPTH];
    MIC_ALIGN mic_f        stack_near_chunk[QUAD_BVH_MAX_STACK_DEPTH];
    MIC_ALIGN unsigned int stack_node_single[QUAD_BVH_MAX_STACK_DEPTH];

    prefetch<PFHINT_L1EX>(stack_node_chunk);      
    prefetch<PFHINT_L1EX>(stack_near_chunk);      

    const mic_m m_active = m_active_i;
    const BVH4AOS* bvh = This->bvh;
    const Node      *const __restrict__ qbvh  = (const Node*) bvh->nodePtr(); //Renderer::scene->qbvh;
    const Triangle1 *const __restrict__ accel = (const Triangle1*) bvh->triPtr(); //Renderer::scene->accel;

    const mic_f min_dist = ray.tnear;
    const mic_f max_dist = ray.tfar;
    
    const mic_f orgx = ray.org.x; //origin.x;
    const mic_f orgy = ray.org.y; //origin.y;
    const mic_f orgz = ray.org.z; //origin.z;  
    
    mic_f max_distance = sel(m_active,max_dist,mic_f::zero());

    const mic_f direction_x = sel(eqz(ray.dir.x),mic_f::ulp(),ray.dir.x);
    const mic_f direction_y = sel(eqz(ray.dir.y),mic_f::ulp(),ray.dir.y);
    const mic_f direction_z = sel(eqz(ray.dir.z),mic_f::ulp(),ray.dir.z);


    const mic_f rdirx = rcp(direction_x);
    const mic_f rdiry = rcp(direction_y);
    const mic_f rdirz = rcp(direction_z);

    const mic_f org_rdirx = orgx * rdirx;
    const mic_f org_rdiry = orgy * rdiry;
    const mic_f org_rdirz = orgz * rdirz;

    const mic_f inf = mic_f::inf();

    /* --  dummy node at index 0, removes branch in inner-loop -- */
    stack_node_chunk[0]  = QBVH_TERMINAL_TOKEN;
    stack_near_chunk[0]  = inf;
  
    unsigned int* __restrict__ sptr_node = stack_node_chunk + 1;
    mic_f       * __restrict__ sptr_near = stack_near_chunk + 1;

    *sptr_node++ = bvh->root; //qbvh[0].lower[0].d;
    *sptr_near++ = sel(m_active,min_dist,inf);

    prefetch<PFHINT_L1EX>(stack_node_single);      

    STAT_COLLECTOR(StatCollector::chunk.numChunks++; StatCollector::chunk.numRays+=countbits(m_active));

    while (1)
      {
	unsigned int curNode = *(sptr_node-1);
	mic_f minDist = upconv16f((const float*)(sptr_near-1));

	mic_m m_dist = lt(minDist,max_distance);

	sptr_node--;
	sptr_near--;

	if (unlikely(curNode == QBVH_TERMINAL_TOKEN)) break;

	STAT_COLLECTOR(StatCollector::chunk.numStackPops++);

	if (unlikely( m_dist == 0)) continue;

	assert(curNode != QBVH_TERMINAL_TOKEN);
	assert(eq(m_dist,minDist,inf) == 0);

	// ==============================
	// === trace single rays here ===
	// ==============================
        if (unlikely(countbits(m_dist) <= SIMD_UTIL_SWITCH_THRESHOLD)) 
        { 
	    float *__restrict__ stack_near_single = (float*)sptr_near;

	    STAT_COLLECTOR(StatCollector::chunk.debug0++);
	    STAT_COLLECTOR(StatCollector::single.numChunks+=countbits(m_dist); StatCollector::single.numRays+=countbits(m_dist));

	    store16f(MIC_M_ALL,&/*ray.tfar*/ray.tfar,max_distance);
	    long rayIndex = -1;
	    const unsigned int m_index = toInt(m_dist);
	    while((rayIndex = bsf64(rayIndex,m_index)) != MIC_NO_BIT_SET_64) 
	      {	    
		stack_node_single[0] = QBVH_TERMINAL_TOKEN;
		stack_node_single[1] = curNode;
		unsigned int sindex = 2;
		// === todo: precompute SOAtoAOS transformation, load with 4x broadcast
		const mic_f org_aos      = SOAtoAOS_4f(rayIndex,ray.org.x,ray.org.y,ray.org.z);
		const mic_f dir_aos      = SOAtoAOS_4f(rayIndex,direction_x,direction_y,direction_z);
		const mic_f rdir_aos     = SOAtoAOS_4f(rayIndex,rdirx,rdiry,rdirz);
		const mic_f org_rdir_aos = org_aos * rdir_aos;
		const mic_f min_dist_aos = upconv1f(&min_dist[rayIndex]); 
      
		//mic_i triangleID         = upconv1i(&ray.triID[rayIndex]);
		mic_f max_dist_aos       = upconv1f(&/*ray.tfar*/ray.tfar[rayIndex]); 
      
		*(mic_f*)stack_near_single = inf; 

		const unsigned int tmp_mask = TMP_QBVH_LEAF_MASK;

		while (1)
		  {
		    unsigned int curNode = stack_node_single[sindex-1];
		    sindex--;

		    STAT_COLLECTOR(StatCollector::single.numStackPops++);
      

		    while (1) 
		      {

			if (unlikely(qbvhLeaf(curNode,tmp_mask))) break;

			const Node* __restrict__ const qptr = qbvhChildPtrNoMask(qbvh,curNode);
			const float* __restrict const b_min = (float*)qptr->lower;
			const float* __restrict const b_max = (float*)qptr->upper;

			prefetch<PFHINT_L1>((char*)qbvh + (unsigned long)curNode);
			prefetch<PFHINT_L1>((char*)qbvh + (unsigned long)curNode + 64);

			mic_f near4_min = inf;
			STAT_COLLECTOR(StatCollector::single.numTraversalSteps++);

			const mic_f clipMinXYZ = upconv16f(b_min) * rdir_aos - org_rdir_aos;
			const mic_f clipMaxXYZ = upconv16f(b_max) * rdir_aos - org_rdir_aos;

			curNode = stack_node_single[sindex-1];

			mic_f clipMinXYZ_min = min_dist_aos;
			clipMinXYZ_min = mask_min(0x7777,clipMinXYZ_min,clipMinXYZ,clipMaxXYZ);
			mic_f clipMaxXYZ_max = max_dist_aos;
			clipMaxXYZ_max = mask_max(0x7777,clipMaxXYZ_max,clipMinXYZ,clipMaxXYZ);

			prefetch<PFHINT_L2>((char*)qbvh + (unsigned long)curNode);
			prefetch<PFHINT_L2>((char*)qbvh + (unsigned long)curNode + 64);

			sindex--;

			const mic_f near4 = set_max4(clipMinXYZ_min);
			const mic_f far4  = set_min4(clipMaxXYZ_max);  

			const mic_m m_lr_hit = le(0x8888,near4,far4);

			store16f(pack_near,near4);


			STAT_COLLECTOR(StatCollector::single.histo[countbits(m_lr_hit)]++);


			near4_min = sel(m_lr_hit,near4,near4_min);

			if (unlikely(ez_mask(m_lr_hit) != 0)) continue;

			const unsigned long i_lr_hit = toInt(m_lr_hit);
			const unsigned long pos_first = bsf64(i_lr_hit);
			const unsigned long num_m_lr_hit = countbits(i_lr_hit); 
	      
			sindex++;
			curNode = ((unsigned int*)b_min)[pos_first];


			if (likely(num_m_lr_hit == 1)) continue;

	      
			const unsigned long pos_sec = bsf64(pos_first,i_lr_hit);


			if (likely(num_m_lr_hit == 2))
			  {
			    const unsigned int dist_first = ((unsigned int*)pack_near)[pos_first];
			    const unsigned int dist_sec   = ((unsigned int*)pack_near)[pos_sec];
			    const unsigned int node_first = curNode;
			    const unsigned int node_sec   = ((unsigned int*)b_min)[pos_sec];

			    if (dist_first <= dist_sec)
			      {
				assert(pack_near[pos_first] <= pack_near[pos_sec]);
				stack_node_single[sindex] = node_sec;
				((unsigned int*)stack_near_single)[sindex] = dist_sec;                      
				sindex++;


				continue;
			      }
			    else
			      {
				assert(pack_near[pos_sec] <= pack_near[pos_first]);
				stack_node_single[sindex] = curNode;
				((unsigned int*)stack_near_single)[sindex] = dist_first;
				curNode = node_sec;
				sindex++;
				continue;
			      }
			  }

			const mic_f child_min_dist = set_min_lanes(near4_min);
			const unsigned int old_sindex = sindex;
			sindex += countbits(i_lr_hit) - 1;

			const mic_m m_child_min_dist = eq(m_lr_hit,child_min_dist,near4);
			const unsigned long pos = bsf64(m_child_min_dist);

			const mic_m m_pos = andn(m_lr_hit,andn(m_child_min_dist,(mic_m)((unsigned int)m_child_min_dist - 1)));
			const mic_i b_min_node = load16i((int*)b_min);

			curNode = ((unsigned int*)b_min)[pos];

			compactustore16i(m_pos,&stack_node_single[old_sindex],b_min_node);
			compactustore16f(m_pos,&stack_near_single[old_sindex],near4);

			// todo: try gather prefetch L2 here !
		      }

		    if (unlikely(curNode == QBVH_TERMINAL_TOKEN)) { break; } // === stack empty ===

		    STAT_COLLECTOR(StatCollector::single.numLeafIntersections++);

		    unsigned int itemOffset    = qbvhItemOffset(curNode);
		    const Triangle1  *__restrict__ tptr    = (Triangle1*)((char*)accel + itemOffset);

		    prefetch<PFHINT_L1>(tptr + 3);
		    prefetch<PFHINT_L1>(tptr + 2);
		    prefetch<PFHINT_L1>(tptr + 1);
		    prefetch<PFHINT_L1>(tptr + 0); // MISS

		    const unsigned int items = qbvhItems(curNode);
		    assert(items > 0);
		    assert(items <= 4);

		    unsigned int triID = qbvhItemOffsetToID(itemOffset);

		    const mic_f zero = mic_f::zero();

		    mic_i triID4 = upconv1i((const int*)&triID);
		    const mic_i and_mask = mic_i::zlc4();

		    STAT_COLLECTOR(StatCollector::single.numAccessedPrimitives+=4);
		    STAT_COLLECTOR(StatCollector::single.numPrimitiveIntersections++);

 
		    const mic_f v0 = gather_4f_zlc(and_mask,
						   (float*)&tptr[0].v0,
						   (float*)&tptr[1].v0,
						   (float*)&tptr[2].v0,
						   (float*)&tptr[3].v0);

		    const mic_f e1 = gather_4f_zlc(and_mask,
						   (float*)&tptr[0].e1,
						   (float*)&tptr[1].e1,
						   (float*)&tptr[2].e1,
						   (float*)&tptr[3].e1);

		    const mic_f e2 = gather_4f_zlc(and_mask,
						   (float*)&tptr[0].e2,
						   (float*)&tptr[1].e2,
						   (float*)&tptr[2].e2,
						   (float*)&tptr[3].e2);
                    const mic_f normal = lcross_zxy(e2,e1);

		    //const mic_f e1 = v1 - v0;
		    //const mic_f e2 = v2 - v0;	     
		    //const mic_f normal = lcross_zxy(e1,e2);
		    const mic_f org = v0 - org_aos;



		    //const mic_f odzxy = msubr231(dir_aos * swizzle(org,_MM_SWIZ_REG_DACB),org,swizzle(dir_aos,_MM_SWIZ_REG_DACB));	      
                    const mic_f odzxy = msubr231(org * swizzle(dir_aos,_MM_SWIZ_REG_DACB), dir_aos, swizzle(org,_MM_SWIZ_REG_DACB));	      
	
		      
		    const mic_f den = ldot3_zxy(dir_aos,normal);	      


		    const mic_f rcp_den = rcp(den);

		    const mic_f uu = ldot3_zxy(e2,odzxy); 
		    const mic_f vv = ldot3_zxy(e1,odzxy); 

		    const mic_f u = uu * rcp_den;
		    const mic_f v = vv * rcp_den;

		    const mic_m valid_u = ge((mic_m)0x1111,u,zero);
		    const mic_m valid_v = ge(valid_u,v,zero);
		    const mic_m m_aperture = le(valid_v,u+v,mic_f::one()); 

		    const mic_f nom = ldot3_zxy(org,normal);

		    if ( unlikely(ez_mask(m_aperture) != 0) ) continue;
		    const mic_f t = rcp_den*nom;
		    const mic_m m_final  = lt(lt(m_aperture,min_dist_aos,t),t,max_dist_aos);

		    max_dist_aos  = sel(m_final,t,max_dist_aos);
		       
		    // ================ Update hit and compact stack ============
		    if (unlikely(nz_mask(m_final)))
		      {
			prefetch<PFHINT_L1EX>(&ray.tfar);  //prefetch<PFHINT_L1EX>(&ray.tfar);      
			prefetch<PFHINT_L1EX>(&ray.id0); //prefetch<PFHINT_L1EX>(&ray.triID);      
			prefetch<PFHINT_L1EX>(&ray.id1); //prefetch<PFHINT_L1EX>(&ray.shaderID);      
			prefetch<PFHINT_L1EX>(&ray.u); //prefetch<PFHINT_L1EX>(&ray.u);      
			prefetch<PFHINT_L1EX>(&ray.v); //prefetch<PFHINT_L1EX>(&ray.v);      
			prefetch<PFHINT_L1EX>(&ray.Ng.x); //prefetch<PFHINT_L1EX>(&ray.gNormal[0]);      
			prefetch<PFHINT_L1EX>(&ray.Ng.y); //prefetch<PFHINT_L1EX>(&ray.gNormal[1]);      
			prefetch<PFHINT_L1EX>(&ray.Ng.z); //prefetch<PFHINT_L1EX>(&ray.gNormal[2]);      

			const mic_f min_dist = set_min16(max_dist_aos);
			const mic_m m_dist = eq(min_dist,max_dist_aos);

			max_dist_aos = min_dist;

			const unsigned long vecIndex = bsf64(m_dist);
			//const unsigned long triIndex = triangleID4[vecIndex];
			const unsigned long triIndex = triID + (vecIndex >> 2);
			const Triangle1  *__restrict__ tri_ptr = &accel[triIndex];

			const mic_m m_final = m_dist^(m_dist & (mic_m)((unsigned long)m_dist - 1));

			assert(countbits(m_final) == 1);
	      
#if  defined(USE_PRECOMPUTED_GNORMAL)
		      const mic_f gnormalz = mic_f(tri_ptr->normal.z);
		      const mic_f gnormalx = mic_f(tri_ptr->normal.x);
		      const mic_f gnormaly = mic_f(tri_ptr->normal.y);

#else	      
		      const mic_f gnormalz = mic_f(normal[vecIndex+0]);
		      const mic_f gnormalx = mic_f(normal[vecIndex+1]);
		      const mic_f gnormaly = mic_f(normal[vecIndex+2]);
#endif

			compactustore16f_low(m_final,&ray.tfar[rayIndex],min_dist); //compactustore16f_low(m_final,&ray.t[rayIndex],min_dist);
			//compactustore16i_low(m_final,&ray.triID[rayIndex],triangleID4);

			compactustore16f_low(m_final,&ray.u[rayIndex],u); //compactustore16f_low(m_final,&ray.u[rayIndex],u4);
			compactustore16f_low(m_final,&ray.v[rayIndex],v); //compactustore16f_low(m_final,&ray.v[rayIndex],v4);

			compactustore16f_low(m_final,&ray.Ng.x[rayIndex],gnormalx); //compactustore16f_low(m_final,&ray.gNormal[0][rayIndex],gnormalx);
			compactustore16f_low(m_final,&ray.Ng.y[rayIndex],gnormaly); //compactustore16f_low(m_final,&ray.gNormal[1][rayIndex],gnormaly);
			compactustore16f_low(m_final,&ray.Ng.z[rayIndex],gnormalz); //compactustore16f_low(m_final,&ray.gNormal[2][rayIndex],gnormalz);

			//ray.shaderID[rayIndex] = tri_ptr->shaderID;
			ray.id0[rayIndex] = tri_ptr->e1.a;
			ray.id1[rayIndex] = tri_ptr->e2.a;

			// === only after max_dist_aos is updated, do a stack compaction ===

			if (likely(sindex >= 2))
			  {
			    //unsigned int old_sindex = sindex;
			    if (likely(sindex < 16))
			      {
				const unsigned int m_num_stack = mic_i::shift1(sindex) - 1;

				const mic_m m_num_stack_low  = toMask(m_num_stack);
				const mic_f snear_low  = upconv16f(stack_near_single + 0);
				const mic_i snode_low  = upconv16i((int*)stack_node_single + 0);
				const mic_m m_stack_compact_low  = le(m_num_stack_low,snear_low,max_dist_aos) | (mic_m)1;
 		  
				compactustore16f_low(m_stack_compact_low,stack_near_single + 0,snear_low);
				compactustore16i_low(m_stack_compact_low,(int*)stack_node_single + 0,snode_low);
				sindex = countbits(m_stack_compact_low);
				assert(sindex < 16);

			      }
			    else if (likely(sindex < 32))
			      {
				const mic_m m_num_stack_high = toMask(mic_i::shift1(sindex-16) - 1); 

				const mic_f snear_low  = upconv16f(stack_near_single + 0);
				const mic_f snear_high = upconv16f(stack_near_single + 16);
				const mic_i snode_low  = upconv16i((int*)stack_node_single + 0);
				const mic_i snode_high = upconv16i((int*)stack_node_single + 16);
				const mic_m m_stack_compact_low  = le(snear_low,max_dist_aos) | (mic_m)1;
				const mic_m m_stack_compact_high = le(m_num_stack_high,snear_high,max_dist_aos);

				compactustore16f(m_stack_compact_low,      stack_near_single + 0,snear_low);
				compactustore16i(m_stack_compact_low,(int*)stack_node_single + 0,snode_low);
				compactustore16f(m_stack_compact_high,      stack_near_single + countbits(m_stack_compact_low),snear_high);
				compactustore16i(m_stack_compact_high,(int*)stack_node_single + countbits(m_stack_compact_low),snode_high);

				sindex = countbits(m_stack_compact_low) + countbits(m_stack_compact_high);

				assert ((unsigned int)m_num_stack_high == (((unsigned int)mic_i::shift1(sindex) - 1) >> 16));
				assert(sindex < 32);
			      }
			    else
			      {
				const mic_m m_num_stack_32 = toMask(mic_i::shift1(sindex-32) - 1); 

				const mic_f snear_0  = upconv16f(stack_near_single + 0);
				const mic_f snear_16 = upconv16f(stack_near_single + 16);
				const mic_f snear_32 = upconv16f(stack_near_single + 32);
				const mic_i snode_0  = upconv16i((int*)stack_node_single + 0);
				const mic_i snode_16 = upconv16i((int*)stack_node_single + 16);
				const mic_i snode_32 = upconv16i((int*)stack_node_single + 32);
				const mic_m m_stack_compact_0  = le(               snear_0 ,max_dist_aos) | (mic_m)1;
				const mic_m m_stack_compact_16 = le(               snear_16,max_dist_aos);
				const mic_m m_stack_compact_32 = le(m_num_stack_32,snear_32,max_dist_aos);

				sindex = 0;
				compactustore16f(m_stack_compact_0,      stack_near_single + sindex,snear_0);
				compactustore16i(m_stack_compact_0,(int*)stack_node_single + sindex,snode_0);
				sindex += countbits(m_stack_compact_0);
				compactustore16f(m_stack_compact_16,      stack_near_single + sindex,snear_16);
				compactustore16i(m_stack_compact_16,(int*)stack_node_single + sindex,snode_16);
				sindex += countbits(m_stack_compact_16);
				compactustore16f(m_stack_compact_32,      stack_near_single + sindex,snear_32);
				compactustore16i(m_stack_compact_32,(int*)stack_node_single + sindex,snode_32);
				sindex += countbits(m_stack_compact_32);

				assert(sindex < 48);

			      }

			  }
		      }
		  }	      

	      }
	    max_distance = ray.tfar/*ray.tfar*/;
	    continue; 
	  }
	// ===========================
	// ===========================
	// ===========================


	assert(curNode != QBVH_TERMINAL_TOKEN);

	const unsigned int tmp_mask = TMP_QBVH_LEAF_MASK;

	while(1)
	  {
	    if (unlikely(qbvhLeaf(curNode,tmp_mask))) break;

	    sptr_node--;
	    sptr_near--;

	    const Node* __restrict__ const qptr = qbvhChildPtrNoMask(qbvh,curNode);

	    STAT_COLLECTOR(StatCollector::chunk.numTraversalSteps++);

	    prefetch<PFHINT_L1>((mic_f*)qptr + 0); 
	    prefetch<PFHINT_L1>((mic_f*)qptr + 1); 

	    STAT_COLLECTOR(StatCollector::chunk.debug2+=countbits(le(m_active,minDist,max_distance)));

	    STAT_COLLECTOR(unsigned int quadhits=0);
	    minDist = *sptr_near;
	    curNode = *sptr_node;


#pragma unroll(4)
	    for (unsigned int i=0;i<4;i++)
	      {
		STAT_COLLECTOR(StatCollector::chunk.numBoxIntersections++);

		const float* __restrict const b_min = (float*)&qptr->lower[i];
		const float* __restrict const b_max = (float*)&qptr->upper[i];

		const mic_f lclipMinX = b_min[0] * rdirx - org_rdirx;
		const mic_f lclipMinY = b_min[1] * rdiry - org_rdiry;
		const mic_f lclipMinZ = b_min[2] * rdirz - org_rdirz;
		const mic_f lclipMaxX = b_max[0] * rdirx - org_rdirx;
		const mic_f lclipMaxY = b_max[1] * rdiry - org_rdiry;
		const mic_f lclipMaxZ = b_max[2] * rdirz - org_rdirz;

		const mic_f lnearP = _max(_max(_min(lclipMinX, lclipMaxX), _min(lclipMinY, lclipMaxY)), _min(lclipMinZ, lclipMaxZ));
		const mic_f lfarP  = _min(_min(_max(lclipMinX, lclipMaxX), _max(lclipMinY, lclipMaxY)), _max(lclipMinZ, lclipMaxZ));      
		const mic_m lhit   = le(_max(lnearP,min_dist), _min(lfarP,max_distance));

		const mic_f boxDist = sel(lhit,lnearP,inf);

		const mic_m m_update_min = lt(boxDist,minDist);

		if (likely(nz_mask(lhit)))
		  {
		    STAT_COLLECTOR(quadhits++);

		    const unsigned int dd = qptr->lower[i].child;

		    sptr_node++;
		    sptr_near++;
		    assert(sptr_node - stack_node_chunk < QUAD_BVH_MAX_STACK_DEPTH);

		    if (nz_mask(m_update_min))
		      {
			*(sptr_node-1) = curNode;
			*(sptr_near-1) = minDist; 

			minDist = boxDist;
			curNode = dd;
		      }
		    else
		      {
			*(sptr_node-1) = dd;
			*(sptr_near-1) = boxDist; 
		      }

		  }	      
	      }
	    STAT_COLLECTOR(StatCollector::chunk.histo[quadhits]++);
	  }

	if (unlikely(curNode == QBVH_TERMINAL_TOKEN)) break;

	STAT_COLLECTOR(StatCollector::chunk.numLeafIntersections++);

	unsigned int itemOffset    = qbvhItemOffset(curNode);
	const Triangle1  *__restrict__ tptr    = (Triangle1*)((char*)accel + itemOffset);
	const unsigned int items = qbvhItems(curNode);
	unsigned int triID = qbvhItemOffsetToID(itemOffset);


	prefetch<PFHINT_L2>((mic_f*)tptr +  0); 
	prefetch<PFHINT_L2>((mic_f*)tptr +  1); 
	prefetch<PFHINT_L2>((mic_f*)tptr +  2); 
	prefetch<PFHINT_L2>((mic_f*)tptr +  3); 

	assert(items > 0);
	assert(items <= 4);

	const mic_f zero = mic_f::zero();
	const mic_m m_active_leaf = ne(minDist,inf);

	STAT_COLLECTOR(StatCollector::chunk.debug1+=countbits(m_active_leaf));

	for(unsigned int i=0;i<items;i++,tptr++,triID++) 

	  {
	    STAT_COLLECTOR(StatCollector::chunk.numAccessedPrimitives++);
	    STAT_COLLECTOR(StatCollector::chunk.numPrimitiveIntersections++);

	    prefetch<PFHINT_L1>((mic_f*)tptr +  1); 

	    const mic_f v0 = upconv4f((float*)&tptr->v0);
            const mic_f e1 = upconv4f((float*)&tptr->e1);
            const mic_f e2 = upconv4f((float*)&tptr->e2);
	    //const mic_f e1 = upconv4f((float*)&tptr->v1) - v0;
	    //const mic_f e2 = upconv4f((float*)&tptr->v2) - v0;

	    //const mic3f org =  mic3f(swAAAA(v0) - orgx,swBBBB(v0) - orgy,swCCCC(v0) - orgz);
	    const mic_f v0_org_x = swAAAA(v0) - orgx;
	    const mic_f v0_org_y = swBBBB(v0) - orgy;
	    const mic_f v0_org_z = swCCCC(v0) - orgz;

#if defined(USE_PRECOMPUTED_GNORMAL)
	    const mic_f nx(tptr->normal.x);
	    const mic_f ny(tptr->normal.y);
	    const mic_f nz(tptr->normal.z);
#else
	    //const mic_f normal = lcross_zxy(e1,e2);
            const mic_f normal = lcross_zxy(e2,e1);
	    const mic_f nz = swAAAA(normal);
	    const mic_f nx = swBBBB(normal);
	    const mic_f ny = swCCCC(normal);
#endif

	    const mic_f nom = madd231(madd231(v0_org_z * nz,v0_org_y, ny),v0_org_x,nx);
	    const mic_f den = madd231(madd231(direction_z * nz,direction_y,ny),direction_x,nx);

	    const mic_m m_sign = lz(den);


	    const mic_f odz = msubr231(direction_x*v0_org_y,v0_org_x,direction_y);
	    const mic_f odx = msubr231(direction_y*v0_org_z,v0_org_y,direction_z);
	    const mic_f ody = msubr231(direction_z*v0_org_x,v0_org_z,direction_x);

	    const mic_f rcp_den = rcp(den);

	    const mic_f uu = zero - madd231(madd231(odz*swCCCC(e2),ody,swBBBB(e2)),odx,swAAAA(e2));
	    const mic_f vv = zero - madd231(madd231(odz*swCCCC(e1),ody,swBBBB(e1)),odx,swAAAA(e1));
	  	 
	    const mic_f u = uu * rcp_den;
	    const mic_f v = vv * rcp_den;

	    const mic_m valid_u = ge(m_active_leaf,u,zero);
	    const mic_m valid_v = ge(valid_u,v,zero);
	    const mic_m m_aperture = le(valid_v,u+v,mic_f::one()); 

	    const mic_f t = rcp_den*nom;

	    if ( unlikely(m_aperture == 0) ) continue;

	    const mic_i i_triID = mic_i(triID);
	    const mic_i ID0 = upconv1i((const int *)&tptr->e1.a);
	    const mic_i ID1 = upconv1i((const int *)&tptr->e2.a);

	    const mic_m m_final  = lt(lt(m_aperture,min_dist,t),t,max_distance);

	    prefetch<PFHINT_L1EX>(&ray.Ng.x); //prefetch<PFHINT_L1EX>(&ray.gNormal[0]);      
	    prefetch<PFHINT_L1EX>(&ray.Ng.y); //prefetch<PFHINT_L1EX>(&ray.gNormal[1]);      
	    prefetch<PFHINT_L1EX>(&ray.Ng.z); //prefetch<PFHINT_L1EX>(&ray.gNormal[2]);      

	    max_distance = sel(m_final,t,max_distance);

	    prefetch<PFHINT_L1EX>(&ray.u); //prefetch<PFHINT_L1EX>(&ray.u);      
	    prefetch<PFHINT_L1EX>(&ray.v); //prefetch<PFHINT_L1EX>(&ray.v);      

	    //const mic_f rden = sel(den_gt_0,rcp_den,-rcp_den);

	    if ( unlikely(m_final == 0) ) continue;	  

	    //todo: less prefetches, use gather for shaderID
	    prefetch<PFHINT_L1EX>(&ray.id0); //prefetch<PFHINT_L1EX>(&ray.triID);      
	    prefetch<PFHINT_L1EX>(&ray.id1); //prefetch<PFHINT_L1EX>(&ray.shaderID);      
	    prefetch<PFHINT_L1EX>(&ray.tfar); //prefetch<PFHINT_L1EX>(&ray.tfar);      

	    store16f(m_final,&ray.Ng.x,nx); //store16f(m_final,&ray.gNormal[0],nx);
	    store16f(m_final,&ray.Ng.y,ny); //store16f(m_final,&ray.gNormal[1],ny);
	    store16f(m_final,&ray.Ng.z,nz); //store16f(m_final,&ray.gNormal[2],nz);

	    store16f(m_final,&ray.u,u); //store16f(m_final,&ray.u,du);
	    store16f(m_final,&ray.v,v); //store16f(m_final,&ray.v,dv);
	    store16i(m_final,&ray.id0,ID0); //store16i(m_final,&ray.shaderID,shdID);
	    store16i(m_final,&ray.id1,ID1); //store16i(m_final,&ray.triID,i_triID);
	  }

	// ---------------------------------------------
	// ---------------------------------------------
      }

    store16f(MIC_M_ALL,&ray.tfar/*ray.tfar*/,max_distance);
  }


  __mmask BVH4AOSTriangle1Intersector16HybridRef::occluded(const BVH4AOSTriangle1Intersector16HybridRef* This, Ray16& ray, const __mmask m_active_i) 
  {
    MIC_ALIGN float        pack_near[16];
    MIC_ALIGN unsigned int stack_node_chunk[QUAD_BVH_MAX_STACK_DEPTH];
    MIC_ALIGN mic_f        stack_near_chunk[QUAD_BVH_MAX_STACK_DEPTH];

    prefetch<PFHINT_L1EX>(stack_node_chunk);      
    prefetch<PFHINT_L1EX>(stack_near_chunk);      

    const mic_m m_active = m_active_i;
    const BVH4AOS* bvh = This->bvh;
    const Node      *const __restrict__ qbvh  = (const Node*) bvh->nodePtr(); //Renderer::scene->qbvh;
    const Triangle1 *const __restrict__ accel = (const Triangle1*) bvh->triPtr(); //Renderer::scene->accel;
  
    //const mic3f &origin    = ray.org;
    //const mic3f &direction = ray.dir;
    const mic_f min_dist = ray.tnear;
    const mic_f max_dist = ray.tfar;

    const mic_f orgx = ray.org.x; //origin.x;
    const mic_f orgy = ray.org.y; //origin.y;
    const mic_f orgz = ray.org.z; //origin.z;  

    mic_f max_distance = sel(m_active,max_dist,mic_f::zero());

    const mic_f direction_x = sel(eqz(ray.dir.x),mic_f::ulp(),ray.dir.x);
    const mic_f direction_y = sel(eqz(ray.dir.y),mic_f::ulp(),ray.dir.y);
    const mic_f direction_z = sel(eqz(ray.dir.z),mic_f::ulp(),ray.dir.z);

    const mic_f rdirx = rcp(direction_x);
    const mic_f rdiry = rcp(direction_y);
    const mic_f rdirz = rcp(direction_z);


    const mic_f org_rdirx = orgx * rdirx;
    const mic_f org_rdiry = orgy * rdiry;
    const mic_f org_rdirz = orgz * rdirz;
    const mic_f inf = mic_f::inf();

    /* --  dummy node at index 0, removes branch in inner-loop -- */

    stack_node_chunk[0]  = QBVH_TERMINAL_TOKEN;
    stack_near_chunk[0]  = inf;

    prefetch<PFHINT_L1EX>(pack_near);      


    unsigned int* __restrict__ sptr_node = stack_node_chunk + 1;
    mic_f       * __restrict__ sptr_near = stack_near_chunk + 1;

    *sptr_node++ = bvh->root; //qbvh[0].lower[0].d;
    *sptr_near++ = sel(m_active,min_dist,inf);

    //

    mic_m m_not_occluded = m_active;

    STAT_COLLECTOR(StatCollector::chunk.numChunks++; StatCollector::chunk.numRays+=countbits(m_active));

    while (1)
      {
	unsigned int curNode = *(sptr_node-1);
	mic_f minDist = upconv16f((const float*)(sptr_near-1));

	mic_m m_dist = lt(m_not_occluded,minDist,max_distance);

	sptr_node--;
	sptr_near--;

	if (unlikely(curNode == QBVH_TERMINAL_TOKEN)) break;

	STAT_COLLECTOR(StatCollector::chunk.numStackPops++);

	if (unlikely( m_dist == 0)) continue;

	assert(curNode != QBVH_TERMINAL_TOKEN);
	assert(eq(m_dist,minDist,inf) == 0);

	// ==============================
	// === trace single rays here ===
	// ==============================
        if (unlikely(countbits(m_dist) <= SIMD_UTIL_SWITCH_THRESHOLD)) 
	  { 
	    unsigned int *__restrict__ stack_node_single = sptr_node;

	    STAT_COLLECTOR(StatCollector::single.numChunks+=countbits(m_dist); StatCollector::single.numRays+=countbits(m_dist));

	    mic_i not_occluded = sel(m_not_occluded,mic_i::minus_one(),mic_i::zero());

	    long rayIndex = -1;

	    const unsigned int m_index = toInt(m_dist);
	    while((rayIndex = bsf64(rayIndex,m_index)) != MIC_NO_BIT_SET_64) 
	      {	    
#if defined(ENABLE_LIGHTSOURCES_MASK)
		const mic_i ray_data = mic_i(ray.data[rayIndex]);
#endif
		stack_node_single[0] = QBVH_TERMINAL_TOKEN;
		stack_node_single[1] = curNode;
		unsigned int sindex = 2;
		// === todo: precompute SOAtoAOS transformation, load with 4x broadcast
		const mic_f org_aos      = SOAtoAOS_4f(rayIndex,ray.org.x,ray.org.y,ray.org.z);
		const mic_f dir_aos      = SOAtoAOS_4f(rayIndex,direction_x,direction_y,direction_z);
		const mic_f rdir_aos     = SOAtoAOS_4f(rayIndex,rdirx,rdiry,rdirz);
		const mic_f org_rdir_aos = org_aos * rdir_aos;
		const mic_f min_dist_aos = upconv1f(&min_dist[rayIndex]);       
		const mic_f max_dist_aos = upconv1f(&max_dist[rayIndex]); 
		const unsigned int tmp_mask = /*TMP_*/QBVH_LEAF_MASK;


		while (1)
		  {
		    unsigned int curNode = stack_node_single[sindex-1];
		    sindex--;

		    STAT_COLLECTOR(StatCollector::single.numStackPops++);
      
		    while (1) 
		      {

			if (unlikely(qbvhLeaf(curNode,tmp_mask))) break;

			const Node* __restrict__ const qptr = qbvhChildPtrNoMask(qbvh,curNode);
			const float* __restrict const b_min = (float*)qptr->lower;
			const float* __restrict const b_max = (float*)qptr->upper;

			prefetch<PFHINT_L1>((char*)qbvh + (unsigned long)curNode);
			prefetch<PFHINT_L1>((char*)qbvh + (unsigned long)curNode + 64);


			mic_f near4_min = inf;
			STAT_COLLECTOR(StatCollector::single.numTraversalSteps++);
			STAT_COLLECTOR(StatCollector::single.numBoxIntersections+=4);

			const mic_f clipMinXYZ = upconv16f(b_min) * rdir_aos - org_rdir_aos;
			const mic_f clipMaxXYZ = upconv16f(b_max) * rdir_aos - org_rdir_aos;

			curNode = stack_node_single[sindex-1];

			mic_f clipMinXYZ_min = min_dist_aos;
			clipMinXYZ_min = mask_min(0x7777,clipMinXYZ_min,clipMinXYZ,clipMaxXYZ);
			mic_f clipMaxXYZ_max = max_dist_aos;
			clipMaxXYZ_max = mask_max(0x7777,clipMaxXYZ_max,clipMinXYZ,clipMaxXYZ);
		      
			prefetch<PFHINT_L2>((char*)qbvh + (unsigned long)curNode); 
			prefetch<PFHINT_L2>((char*)qbvh + (unsigned long)curNode + 64);

			sindex--;

			const mic_f near4 = set_max4(clipMinXYZ_min);
			const mic_f far4  = set_min4(clipMaxXYZ_max);  

			const mic_m m_lr_hit = le(0x8888,near4,far4);

			store16f(pack_near,near4);


			STAT_COLLECTOR(StatCollector::single.histo[countbits(m_lr_hit)]++);


			near4_min = sel(m_lr_hit,near4,near4_min);

			if (unlikely(ez_mask(m_lr_hit) != 0)) continue;

			const unsigned long i_lr_hit = toInt(m_lr_hit);
			const unsigned long pos_first = bsf64(i_lr_hit);
			const unsigned long num_m_lr_hit = countbits(i_lr_hit); 
	      
			sindex++;
			curNode = ((unsigned int*)b_min)[pos_first];


			if (likely(num_m_lr_hit == 1)) continue;

	      
			const unsigned long pos_sec = bsf64(pos_first,i_lr_hit);


			if (likely(num_m_lr_hit == 2))
			  {
			    const unsigned int dist_first = ((unsigned int*)pack_near)[pos_first];
			    const unsigned int dist_sec   = ((unsigned int*)pack_near)[pos_sec];
			    const unsigned int node_first = curNode;
			    const unsigned int node_sec   = ((unsigned int*)b_min)[pos_sec];

			    if (dist_first <= dist_sec)
			      {
				assert(pack_near[pos_first] <= pack_near[pos_sec]);
				stack_node_single[sindex] = node_sec;
				sindex++;


				continue;
			      }
			    else
			      {
				assert(pack_near[pos_sec] <= pack_near[pos_first]);
				stack_node_single[sindex] = curNode;
				curNode = node_sec;
				sindex++;
				continue;
			      }
			  }

			const mic_f child_min_dist = set_min_lanes(near4_min);
			const unsigned int old_sindex = sindex;
			sindex += countbits(i_lr_hit) - 1;

			const mic_m m_child_min_dist = eq(m_lr_hit,child_min_dist,near4);
			const unsigned long pos = bsf64(m_child_min_dist);

			const mic_m m_pos = andn(m_lr_hit,andn(m_child_min_dist,(mic_m)((unsigned int)m_child_min_dist - 1)));
			const mic_i b_min_node = load16i((int*)b_min);

			curNode = ((unsigned int*)b_min)[pos];

			compactustore16i(m_pos,&stack_node_single[old_sindex],b_min_node);
			// todo: try gather prefetch L2 here !
		      }

		    if (unlikely(curNode == QBVH_TERMINAL_TOKEN)) { break; } // === stack empty ===

		    STAT_COLLECTOR(StatCollector::single.numLeafIntersections++);

		    unsigned int itemOffset    = qbvhItemOffset(curNode);
		    const Triangle1  *__restrict__ tptr    = (Triangle1*)((char*)accel + itemOffset);

		    prefetch<PFHINT_NT>(tptr + 3);
		    prefetch<PFHINT_NT>(tptr + 2);
		    prefetch<PFHINT_NT>(tptr + 1);
		    prefetch<PFHINT_NT>(tptr + 0); // MISS



		    const unsigned int items = qbvhItems(curNode);
		    assert(items > 0);
		    assert(items <= 4);

		    const mic_f zero = mic_f::zero();
		    const mic_i and_mask = mic_i::zlc4();

		    {		      
		      STAT_COLLECTOR(StatCollector::single.numAccessedPrimitives+=4);
		      STAT_COLLECTOR(StatCollector::single.numPrimitiveIntersections++);

		      const mic_f v0 = gather_4f_zlc(and_mask,
						     (float*)&tptr[0].v0,
						     (float*)&tptr[1].v0,
						     (float*)&tptr[2].v0,
						     (float*)&tptr[3].v0);

		      const mic_f e1 = gather_4f_zlc(and_mask,
						     (float*)&tptr[0].e1,
						     (float*)&tptr[1].e1,
						     (float*)&tptr[2].e1,
						     (float*)&tptr[3].e1);

		      const mic_f e2 = gather_4f_zlc(and_mask,
						     (float*)&tptr[0].e2,
						     (float*)&tptr[1].e2,
						     (float*)&tptr[2].e2,
						     (float*)&tptr[3].e2);

                      const mic_f normal = lcross_zxy(e2,e1);
		      //const mic_f e1 = v1 - v0;
		      //const mic_f e2 = v2 - v0;	     
		      //const mic_f normal = lcross_zxy(e1,e2);
		      const mic_f org = v0 - org_aos;


		      //const mic_f odzxy = msubr231(dir_aos * swizzle(org,_MM_SWIZ_REG_DACB),org,swizzle(dir_aos,_MM_SWIZ_REG_DACB));
                      const mic_f odzxy = msubr231(org * swizzle(dir_aos,_MM_SWIZ_REG_DACB), dir_aos, swizzle(org,_MM_SWIZ_REG_DACB));	      
                      const mic_f den = ldot3_zxy(dir_aos,normal);
 
		      const mic_f uu = ldot3_zxy(e2,odzxy); 
		      const mic_f rcp_den = rcp(den);

		      const mic_f vv = ldot3_zxy(e1,odzxy); 


		      const mic_f u = uu * rcp_den;
		      const mic_f v = vv * rcp_den;

		      const mic_m valid_u = ge((mic_m)0x1111,u,zero);
		      const mic_m valid_v = ge(valid_u,v,zero);
		      const mic_m m_aperture = le(valid_v,u+v,mic_f::one()); 

		      const mic_f nom = ldot3_zxy(org,normal);

		      const mic_f t = rcp_den*nom;

		      if ( unlikely(ez_mask(m_aperture) != 0) ) continue;

		      mic_m m_final  = lt(lt(m_aperture,min_dist_aos,t),t,max_dist_aos);

#if defined(ENABLE_LIGHTSOURCES_MASK)
		      const mic_i lsMask0 = mic_i(tptr[0].data[0]);
		      const mic_i lsMask1 = mic_i(tptr[1].data[0]);
		      const mic_i lsMask2 = mic_i(tptr[2].data[0]);
		      const mic_i lsMask3 = mic_i(tptr[3].data[0]);
		      mic_i lsMask = lsMask0;
		      lsMask = sel(0xf0  ,lsMask1,lsMask);
		      lsMask = sel(0xf00 ,lsMask2,lsMask);
		      lsMask = sel(0xf000,lsMask3,lsMask);
		      m_final = mask_test(m_final,lsMask,ray_data);

#endif	      
		      if (unlikely(m_final != 0))
			{
			  not_occluded[rayIndex] = 0;
			  break;
			}

		      //m_update |= m_final;
		    }
	  
		  }

	      }
	    m_not_occluded = eq(not_occluded,mic_i::minus_one());
	    if (unlikely(m_not_occluded == 0)) {
	      return m_active;
	    }
	    continue; 
	  }
	// ===========================
	// ===========================
	// ===========================


	assert(curNode != QBVH_TERMINAL_TOKEN);

	const unsigned int tmp_mask = /*TMP_*/QBVH_LEAF_MASK;

	while(1)
	  {
	    if (unlikely(qbvhLeaf(curNode,tmp_mask))) break;

	    sptr_node--;
	    sptr_near--;

	    const Node* __restrict__ const qptr = qbvhChildPtrNoMask(qbvh,curNode);

	    STAT_COLLECTOR(StatCollector::chunk.numTraversalSteps++);
		
	    prefetch<PFHINT_L1>((mic_f*)qptr + 0); 
	    prefetch<PFHINT_L1>((mic_f*)qptr + 1); 

	    STAT_COLLECTOR(StatCollector::chunk.debug2+=countbits(le(m_active,minDist,max_distance)));

	    STAT_COLLECTOR(unsigned int quadhits=0);
	    minDist = *sptr_near;
	    curNode = *sptr_node;


#pragma unroll(4)
	    for (unsigned int i=0;i<4;i++)
	      {
		STAT_COLLECTOR(StatCollector::chunk.numBoxIntersections++);

		const float* __restrict const b_min = (float*)&qptr->lower[i];
		const float* __restrict const b_max = (float*)&qptr->upper[i];

		const mic_f lclipMinX = b_min[0] * rdirx - org_rdirx;
		const mic_f lclipMinY = b_min[1] * rdiry - org_rdiry;
		const mic_f lclipMinZ = b_min[2] * rdirz - org_rdirz;
		const mic_f lclipMaxX = b_max[0] * rdirx - org_rdirx;
		const mic_f lclipMaxY = b_max[1] * rdiry - org_rdiry;
		const mic_f lclipMaxZ = b_max[2] * rdirz - org_rdirz;

		const mic_f lnearP = _max(_max(_min(lclipMinX, lclipMaxX), _min(lclipMinY, lclipMaxY)), _min(lclipMinZ, lclipMaxZ));
		const mic_f lfarP  = _min(_min(_max(lclipMinX, lclipMaxX), _max(lclipMinY, lclipMaxY)), _max(lclipMinZ, lclipMaxZ));      
		const mic_m lhit   = le(m_not_occluded,_max(lnearP,min_dist), _min(lfarP,max_distance));

		const mic_f boxDist = sel(lhit,lnearP,inf);

		const mic_m m_update_min = lt(m_not_occluded,boxDist,minDist);

		if (likely(nz_mask(lhit)))
		  {
		    STAT_COLLECTOR(quadhits++);

		    const unsigned int dd = qptr->lower[i].child;

		    sptr_node++;
		    sptr_near++;
		    assert(sptr_node - stack_node_chunk < QUAD_BVH_MAX_STACK_DEPTH);

		    if (nz_mask(m_update_min))
		      {
			*(sptr_node-1) = curNode;
			*(sptr_near-1) = minDist; 

			minDist = boxDist;
			curNode = dd;
		      }
		    else
		      {

			*(sptr_node-1) = dd;
			*(sptr_near-1) = boxDist; 
		      }

		  }	      
	      }

	    STAT_COLLECTOR(StatCollector::chunk.histo[quadhits]++);
	  }

	if (unlikely(curNode == QBVH_TERMINAL_TOKEN)) break;    

	STAT_COLLECTOR(StatCollector::chunk.numLeafIntersections++);

	unsigned int itemOffset    = qbvhItemOffset(curNode);
	const Triangle1  *__restrict__ tptr    = (Triangle1*)((char*)accel + itemOffset);
	const unsigned int items = qbvhItems(curNode);

	prefetch<PFHINT_NT>((mic_f*)tptr +  0); 
	prefetch<PFHINT_L2>((mic_f*)tptr +  1); 
	prefetch<PFHINT_L2>((mic_f*)tptr +  2); 
	prefetch<PFHINT_L2>((mic_f*)tptr +  3); 

	assert(items > 0);
	assert(items <= 4);

	const mic_f zero = mic_f::zero();
	mic_m m_active_leaf = ne(m_not_occluded,minDist,inf);

	STAT_COLLECTOR(StatCollector::chunk.debug1+=countbits(m_active_leaf));

	for(unsigned int i=0;i<items;i++,tptr++) 
	  {
	    STAT_COLLECTOR(StatCollector::chunk.numAccessedPrimitives++);
	    STAT_COLLECTOR(StatCollector::chunk.numPrimitiveIntersections++);

	    prefetch<PFHINT_NT>((mic_f*)tptr +  1); 

	    const mic_f v0 = upconv4f((float*)&tptr->v0);
            const mic_f e1 = upconv4f((float*)&tptr->e1);
            const mic_f e2 = upconv4f((float*)&tptr->e2);
	    //const mic_f e1 = upconv4f((float*)&tptr->v1) - v0;
	    //const mic_f e2 = upconv4f((float*)&tptr->v2) - v0;

	    //const mic3f org =  mic3f(swAAAA(v0) - orgx,swBBBB(v0) - orgy,swCCCC(v0) - orgz);
	    const mic_f v0_org_x = swAAAA(v0) - orgx;
	    const mic_f v0_org_y = swBBBB(v0) - orgy;
	    const mic_f v0_org_z = swCCCC(v0) - orgz;

#if defined(USE_PRECOMPUTED_GNORMAL)
	    const mic_f nx(tptr->normal.x);
	    const mic_f ny(tptr->normal.y);
	    const mic_f nz(tptr->normal.z);
#else
	    //const mic_f normal = lcross_zxy(e1,e2);
            const mic_f normal = lcross_zxy(e2,e1);
	    const mic_f nz = swAAAA(normal);
	    const mic_f nx = swBBBB(normal);
	    const mic_f ny = swCCCC(normal);
#endif

	    const mic_f den = madd231(madd231(direction_z * nz,direction_y,ny),direction_x,nx);

	    const mic_f odz = msubr231(direction_x*v0_org_y,v0_org_x,direction_y);
	    const mic_f odx = msubr231(direction_y*v0_org_z,v0_org_y,direction_z);
	    const mic_f ody = msubr231(direction_z*v0_org_x,v0_org_z,direction_x);

	    const mic_f rcp_den = rcp(den);

	    const mic_f uu = zero - madd231(madd231(odz*swCCCC(e2),ody,swBBBB(e2)),odx,swAAAA(e2));
	    const mic_f vv = zero - madd231(madd231(odz*swCCCC(e1),ody,swBBBB(e1)),odx,swAAAA(e1));
	  	 
	    const mic_f u = uu * rcp_den;
	    const mic_f v = vv * rcp_den;

	    const mic_m valid_u = ge(m_active_leaf,u,zero);
	    const mic_m valid_v = ge(valid_u,v,zero);
	    const mic_m m_aperture = le(valid_v,u+v,mic_f::one()); 


	    const mic_f nom = madd231(madd231(v0_org_z * nz,v0_org_y, ny),v0_org_x,nx);
	  
	    if ( unlikely(m_aperture == 0) ) continue;

	    const mic_f t = rcp_den*nom;

	    mic_m m_final  = lt(lt(m_aperture,min_dist,t),t,max_distance);

#if defined(ENABLE_LIGHTSOURCES_MASK)
	    m_final = mask_test(m_final,ray.data,upconv1i((int*)&tptr->data[0]));

#endif	      

	    max_distance = sel(m_final,t,max_distance);

	    m_active_leaf  &= ~m_final;
	    m_not_occluded &= ~m_final;
	    if (unlikely(m_active_leaf == 0)) break;
	    //if ( unlikely(m_final == 0) ) continue;

	  }
	if (unlikely(m_not_occluded == 0)) {
	  return m_active;
	}

	// ---------------------------------------------
	// ---------------------------------------------
      }
    return m_active & ~m_not_occluded;
  }

  void BVH4AOSTriangle1Intersector16HybridRefRegister () {
    TriangleMesh::intersectors16.add("bvh4aos","triangle1","hybrid_ref","moeller",true ,BVH4AOSTriangle1Intersector16HybridRef::create);
    TriangleMesh::intersectors16.setAccelDefaultTraverser("bvh4aos","hybrid_ref");
  }
}
