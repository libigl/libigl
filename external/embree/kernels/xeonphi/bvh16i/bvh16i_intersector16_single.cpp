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

#include "bvh16i_intersector16_single.h"
#include "geometry/triangle1.h"
#include "bvh16i.h"


//#define NEAR_FAR_OPT_INTERSECT
//#define NEAR_FAR_OPT_OCCLUDED

#define ENABLE_EXTENDED_LEAVES

__forceinline void sort(void * stack_dist,
			void * stack_node,
			const size_t sindex0,
			const size_t sindex1)
{
  if (unlikely(((unsigned int*)stack_dist)[sindex1] < ((unsigned int*)stack_dist)[sindex0]))
    {
      std::swap(((unsigned int*)stack_dist)[sindex1],((unsigned int*)stack_dist)[sindex0]);
      std::swap(((unsigned int*)stack_node)[sindex1],((unsigned int*)stack_node)[sindex0]);
    }
 
}

		  // const BVH16i::Node* __restrict__ const next = (BVH16i::Node*)((char*)bvh16 + curNode.id());

		  // prefetch<PFHINT_L2>(&next->min_x);
		  // prefetch<PFHINT_L2>(&next->max_x);
		  // prefetch<PFHINT_L2>(&next->min_y);
		  // prefetch<PFHINT_L2>(&next->max_y);
		  // prefetch<PFHINT_L2>(&next->min_z);
		  // prefetch<PFHINT_L2>(&next->max_z);
		  // prefetch<PFHINT_L2>(&next->child);
			
namespace embree
{
  namespace isa
  {

    static unsigned int TMP_BVH16I_LEAF_MASK = BVH16I_LEAF_MASK; // needed due to compiler efficiency bug

    static __aligned(64) int zlc4[4] = {0xffffffff,0xffffffff,0xffffffff,0};

    void BVH16iIntersector16Single::intersect(mic_i* valid_i, BVH16i* bvh, Ray16& ray16)
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

      const Triangle1 * __restrict__ accel = (Triangle1*)bvh->triPtr();
      BVH16i::Node* __restrict__ bvh16 = bvh->node16;

      stack_node[0] = BVH16I_TERMINAL_TOKEN;
      
      NodeRef root = bvh16->child[0];

      const unsigned int leaf_mask = TMP_BVH16I_LEAF_MASK;

#if defined(NEAR_FAR_OPT_INTERSECT)            
      const mic_i nearXoffset16 = select(rdir16.x >= 0.0f,mic_i(0),mic_i(sizeof(mic_f)));
      const mic_i nearYoffset16 = select(rdir16.y >= 0.0f,mic_i(0),mic_i(sizeof(mic_f)));
      const mic_i nearZoffset16 = select(rdir16.z >= 0.0f,mic_i(0),mic_i(sizeof(mic_f)));
#endif

      long rayIndex = -1;
      while((rayIndex = bitscan64(rayIndex,toInt(m_valid))) != BITSCAN_NO_BIT_SET_64)	    
        {
	  stack_node[1] = root;
	  size_t sindex = 2;

	  const mic_f dir_xyz  = loadAOS4to16f(rayIndex,ray16.dir.x,ray16.dir.y,ray16.dir.z);
	  const mic_f rdir_xyz = rcp(dir_xyz);
	  const mic_f org_xyz  = loadAOS4to16f(rayIndex,ray16.org.x,ray16.org.y,ray16.org.z);
	  const mic_f org_rdir_xyz = rdir_xyz * org_xyz;
	  const mic3f rdir(rdir16.x[rayIndex],rdir16.y[rayIndex],rdir16.z[rayIndex]);
	  const mic3f org(swAAAA(org_xyz),swBBBB(org_xyz),swCCCC(org_xyz));
	  const mic3f org_rdir = rdir * org;	  

	  //const mic3f org_rdir(swAAAA(org_rdir_xyz),swBBBB(org_rdir_xyz),swCCCC(org_rdir_xyz));

	  const mic_f tnear = broadcast1to16f(&ray16.tnear[rayIndex]);
	  mic_f       tfar  = broadcast1to16f(&ray16.tfar[rayIndex]);
	  
	  while (1)
	    {
	      NodeRef curNode = stack_node[sindex-1];
	      sindex--;

#if defined(NEAR_FAR_OPT_INTERSECT)            
	      const unsigned int nearXoffset0 = nearXoffset16[rayIndex];
	      const unsigned int nearYoffset0 = nearYoffset16[rayIndex];
	      const unsigned int nearZoffset0 = nearZoffset16[rayIndex];

	      const unsigned int nearXoffset1 = nearXoffset0 ^ sizeof(mic_f);
	      const unsigned int nearYoffset1 = nearYoffset0 ^ sizeof(mic_f);
	      const unsigned int nearZoffset1 = nearZoffset0 ^ sizeof(mic_f);
#endif
	      while (1) 
		{
		  /* test if this is a leaf node */
		  if (unlikely(curNode.isLeaf(leaf_mask))) break;       
		  STAT3(normal.trav_nodes,1,1,1);

		  const BVH16i::Node* __restrict__ const bptr = (BVH16i::Node*)((char*)bvh16 + curNode.id());

		  prefetch<PFHINT_L1>(&bptr->min_x);
		  prefetch<PFHINT_L1>(&bptr->max_x);
		  prefetch<PFHINT_L1>(&bptr->min_y);
		  prefetch<PFHINT_L1>(&bptr->max_y);
		  prefetch<PFHINT_L1>(&bptr->min_z);
		  prefetch<PFHINT_L1>(&bptr->max_z);
		  prefetch<PFHINT_L1>(&bptr->child);


#if !defined(NEAR_FAR_OPT_INTERSECT)
		  const mic_f min_x = bptr->min_x * rdir.x - org_rdir.x;
		  const mic_f max_x = bptr->max_x * rdir.x - org_rdir.x;

		  const mic_f min_y = bptr->min_y * rdir.y - org_rdir.y;
		  const mic_f max_y = bptr->max_y * rdir.y - org_rdir.y;

		  const mic_f min_z = bptr->min_z * rdir.z - org_rdir.z;
		  const mic_f max_z = bptr->max_z * rdir.z - org_rdir.z;

		  const mic_f nearX = min(min_x,max_x);
		  const mic_f farX  = max(min_x,max_x);

		  const mic_f nearY = min(min_y,max_y);
		  const mic_f farY  = max(min_y,max_y);

		  const mic_f nearZ = min(min_z,max_z);
		  const mic_f farZ  = max(min_z,max_z);

#else


		  const mic_f nearX = load16f((float*)((const char*)&bptr->min_x + (size_t)nearXoffset0)) * rdir.x - org_rdir.x;
		  const mic_f farX  = load16f((float*)((const char*)&bptr->min_x + (size_t)nearXoffset1)) * rdir.x - org_rdir.x;

		  const mic_f nearY = load16f((float*)((const char*)&bptr->min_y + (size_t)nearYoffset0)) * rdir.y - org_rdir.y;
		  const mic_f farY  = load16f((float*)((const char*)&bptr->min_y + (size_t)nearYoffset1)) * rdir.y - org_rdir.y;

		  const mic_f nearZ = load16f((float*)((const char*)&bptr->min_z + (size_t)nearZoffset0)) * rdir.z - org_rdir.z;
		  const mic_f farZ  = load16f((float*)((const char*)&bptr->min_z + (size_t)nearZoffset1)) * rdir.z - org_rdir.z;
		  
#endif
		  sindex--;
		  curNode = stack_node[sindex]; // early pop of next node

		  const mic_f near16 = max(max(nearX,nearY),max(nearZ,tnear));
		  const mic_f far16  = min(min(farX,farY),min(farZ,tfar));


		  const mic_m hitm = le(near16,far16);
		  const mic_f tNear_pos = select(hitm,near16,inf);

		  STAT3(normal.trav_hit_boxes[countbits(hitm)],1,1,1);

		  /* if no child is hit, continue with early popped child */
		  if (unlikely(none(hitm))) continue;

		  sindex++;
        
		  const unsigned long hiti = toInt(hitm);
		  const unsigned long pos_first = bitscan64(hiti);
		  const unsigned long num_hitm = countbits(hiti); 
        
		  /* if a single child is hit, continue with that child */
		  curNode = bptr->child[pos_first];
		  if (likely(num_hitm == 1)) continue;
        
		  /* if two children are hit, push in correct order */
		  const unsigned long pos_second = bitscan64(pos_first,hiti);
		  if (likely(num_hitm == 2))
		    {
		      const unsigned int dist_first  = ((unsigned int*)&near16)[pos_first];
		      const unsigned int dist_second = ((unsigned int*)&near16)[pos_second];
		      const unsigned int node_first  = curNode;
		      const unsigned int node_second = bptr->child[pos_second];
          
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
#if 0
		  else if (likely(num_hitm <= 4))
		    {
		      const unsigned int old_sindex = sindex;
		      compactustore16f(hitm,&stack_dist[old_sindex],near16);
		      const mic_i children = bptr->child;
		      compactustore16i(hitm,&stack_node[old_sindex],children);
		      sindex += num_hitm;

		      sort(stack_dist,stack_node,sindex-4,sindex-3);
		      sort(stack_dist,stack_node,sindex-2,sindex-1);
		      sort(stack_dist,stack_node,sindex-3,sindex-1);
		      curNode = stack_node[sindex-1];
		      sindex--;
		      continue;
		    }
#endif



		  /* continue with closest child and push all others */
		  const mic_f min_dist = set_min16(tNear_pos);

		  const unsigned int old_sindex = sindex;
		  sindex += countbits(hiti) - 1;
		  assert(sindex < 3*BVH4i::maxDepth+1);
		  const mic_i children = bptr->child;
        
		  const mic_m closest_child = eq(hitm,min_dist,near16);
		  const unsigned long closest_child_pos = bitscan64(closest_child);
		  const mic_m m_pos = andn(hitm,andn(closest_child,(mic_m)((unsigned int)closest_child - 1)));
		  curNode = bptr->child[closest_child_pos];
		  compactustore16f(m_pos,&stack_dist[old_sindex],near16);
		  compactustore16i(m_pos,&stack_node[old_sindex],children);

		  if (unlikely(((unsigned int*)stack_dist)[sindex-2] < ((unsigned int*)stack_dist)[sindex-1]))
		    {
		      std::swap(((unsigned int*)stack_dist)[sindex-2],((unsigned int*)stack_dist)[sindex-1]);
		      std::swap(((unsigned int*)stack_node)[sindex-2],((unsigned int*)stack_node)[sindex-1]);
		    }

		}

	      /* return if stack is empty */
	      if (unlikely(curNode == BVH16I_TERMINAL_TOKEN)) break;
	  
#if defined(ENABLE_EXTENDED_LEAVES)
	      if (likely((curNode.isLeaf(BVH16I_EXTENDED_LEAF_MASK) == BVH16I_EXTENDED_LEAF_MASK)))
		{
		  curNode.id() ^= BVH16I_EXTENDED_LEAF_MASK;

		  const BVH4i::Node* __restrict__ const node4 = (BVH4i::Node*)((char*)bvh16 + curNode.id());

		  const float* __restrict const plower = (float*)node4->lower;
		  const float* __restrict const pupper = (float*)node4->upper;

		  prefetch<PFHINT_L1>((char*)node4 + 0);
		  prefetch<PFHINT_L1>((char*)node4 + 64);
        
		  /* intersect single ray with 4 bounding boxes */
		  const mic_f tLowerXYZ = load16f(plower) * rdir_xyz - org_rdir_xyz;
		  const mic_f tUpperXYZ = load16f(pupper) * rdir_xyz - org_rdir_xyz;
		  const mic_f tLower = mask_min(0x7777,tnear,tLowerXYZ,tUpperXYZ);
		  const mic_f tUpper = mask_max(0x7777,tfar,tLowerXYZ,tUpperXYZ);

		  const mic_f tNear = vreduce_max4(tLower);
		  const mic_f tFar  = vreduce_min4(tUpper);  
		  const mic_m hitm = le(0x8888,tNear,tFar);
		  const mic_f tNear_pos = select(hitm,tNear,inf);

		  STAT3(normal.trav_hit_boxes[countbits(hitm)],1,1,1);

		  /* if no child is hit, continue with early popped child */
		  if (unlikely(none(hitm))) continue;
        
		  const unsigned long hiti = toInt(hitm);
		  const unsigned long pos_first = bitscan64(hiti);
		  const unsigned long num_hitm = countbits(hiti); 
        
		  /* if a single child is hit, continue with that child */
		  curNode = ((unsigned int *)plower)[pos_first];
		  if (likely(num_hitm == 1)) 
		    {
		      
		    }
		  else if (likely(num_hitm == 2)) /* if two children are hit, push in correct order */
		    {
		      const unsigned long pos_second = bitscan64(pos_first,hiti);
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
			}
		      else
			{
			  stack_node[sindex] = curNode;
			  ((unsigned int*)stack_dist)[sindex] = dist_first;
			  curNode = node_second;
			  sindex++;
			  assert(sindex < 3*BVH4i::maxDepth+1);
			}
		    }
		  else
		    {        
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
#endif


	      STAT3(normal.trav_leaves,1,1,1);
	      STAT3(normal.trav_prims,4,4,4);


	      /* intersect one ray against four triangles */

	      //////////////////////////////////////////////////////////////////////////////////////////////////

	      const Triangle1* tptr  = (Triangle1*) curNode.leaf(accel);

	      
	      prefetch<PFHINT_L1>(tptr + 3);
	      prefetch<PFHINT_L1>(tptr + 2);
	      prefetch<PFHINT_L1>(tptr + 1);
	      prefetch<PFHINT_L1>(tptr + 0); 

	      const mic_i and_mask = broadcast4to16i(zlc4);
	      
	      const mic_f v0 = gather_4f_zlc(and_mask,
					     (float*)&tptr[0].v0,
					     (float*)&tptr[1].v0,
					     (float*)&tptr[2].v0,
					     (float*)&tptr[3].v0);
	      
	      const mic_f v1 = gather_4f_zlc(and_mask,
					     (float*)&tptr[0].v1,
					     (float*)&tptr[1].v1,
					     (float*)&tptr[2].v1,
					     (float*)&tptr[3].v1);
	      
	      const mic_f v2 = gather_4f_zlc(and_mask,
					     (float*)&tptr[0].v2,
					     (float*)&tptr[1].v2,
					     (float*)&tptr[2].v2,
					     (float*)&tptr[3].v2);
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

	      const mic_m m_final  = lt(lt(m_aperture,tnear,t),t,tfar);

	      tfar  = select(m_final,t,tfar); 
		    
	      //////////////////////////////////////////////////////////////////////////////////////////////////


	      /* did the ray hot one of the four triangles? */
	      if (unlikely(any(m_final)))
		{
		  STAT3(normal.trav_prim_hits,1,1,1);

		  const mic_f min_dist_new = vreduce_min(tfar);
		  const mic_m m_dist = eq(min_dist_new,tfar);

		  const size_t vecIndex = bitscan(toInt(m_dist));
		  const size_t triIndex = vecIndex >> 2;

		  const Triangle1  *__restrict__ tri_ptr = tptr + triIndex;

		  const mic_m m_tri = m_dist^(m_dist & (mic_m)((unsigned int)m_dist - 1));

		  const mic_f gnormalx = mic_f(tri_ptr->Ng.x);
		  const mic_f gnormaly = mic_f(tri_ptr->Ng.y);
		  const mic_f gnormalz = mic_f(tri_ptr->Ng.z);

#if defined(__USE_RAY_MASK__)
		  if ( (tri_ptr->mask() & ray16.mask[rayIndex]) != 0 )
#else
		  if (1)
#endif
		    {
		      prefetch<PFHINT_L1EX>(&ray16.tfar);  
		      prefetch<PFHINT_L1EX>(&ray16.u);
		      prefetch<PFHINT_L1EX>(&ray16.v);
		      prefetch<PFHINT_L1EX>(&ray16.Ng.x); 
		      prefetch<PFHINT_L1EX>(&ray16.Ng.y); 
		      prefetch<PFHINT_L1EX>(&ray16.Ng.z); 
		      prefetch<PFHINT_L1EX>(&ray16.geomID);
		      prefetch<PFHINT_L1EX>(&ray16.primID);

		      tfar = min_dist_new;
		  
		      compactustore16f_low(m_tri,&ray16.tfar[rayIndex],tfar);
		      compactustore16f_low(m_tri,&ray16.u[rayIndex],u); 
		      compactustore16f_low(m_tri,&ray16.v[rayIndex],v); 
		      compactustore16f_low(m_tri,&ray16.Ng.x[rayIndex],gnormalx); 
		      compactustore16f_low(m_tri,&ray16.Ng.y[rayIndex],gnormaly); 
		      compactustore16f_low(m_tri,&ray16.Ng.z[rayIndex],gnormalz); 

		      ray16.geomID[rayIndex] = tri_ptr->geomID();
		      ray16.primID[rayIndex] = tri_ptr->primID();

		      /* compact the stack if size of stack >= 2 */
		      if (likely(sindex >= 2))
			{
			  if (likely(sindex < 16))
			    {
			      const unsigned int m_num_stack = mic_m::shift1[sindex] - 1;
			      const mic_m m_num_stack_low  = toMask(m_num_stack);
			      const mic_f snear_low  = load16f(stack_dist + 0);
			      const mic_i snode_low  = load16i((int*)stack_node + 0);
			      const mic_m m_stack_compact_low  = le(m_num_stack_low,snear_low,tfar) | (mic_m)1;
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
			      const mic_m m_stack_compact_low  = le(snear_low,tfar) | (mic_m)1;
			      const mic_m m_stack_compact_high = le(m_num_stack_high,snear_high,tfar);
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
			      const mic_m m_stack_compact_0  = le(               snear_0 ,tfar) | (mic_m)1;
			      const mic_m m_stack_compact_16 = le(               snear_16,tfar);
			      const mic_m m_stack_compact_32 = le(m_num_stack_32,snear_32,tfar);

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
		}
	    }	  
	}
    }
    
    void BVH16iIntersector16Single::occluded(mic_i* valid_i, BVH16i* bvh, Ray16& ray16)
    {
      /* near and node stack */
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

      /* setup */
      const mic_m m_valid     = *(mic_i*)valid_i != mic_i(0);
      const mic3f rdir16      = rcp_safe(ray16.dir);
      mic_m terminated        = !m_valid;
      const mic_f inf         = mic_f(pos_inf);
      const mic_f zero        = mic_f::zero();

      const Triangle1 * __restrict__ accel = (Triangle1*)bvh->triPtr();
      BVH16i::Node* __restrict__ bvh16 = bvh->node16;

      stack_node[0] = BVH16I_TERMINAL_TOKEN;
      NodeRef root = bvh16->child[0];

#if defined(NEAR_FAR_OPT_OCCLUDED)            
      const mic_i nearXoffset16 = select(rdir16.x >= 0.0f,mic_i(0),mic_i(sizeof(mic_f)));
      const mic_i nearYoffset16 = select(rdir16.y >= 0.0f,mic_i(0),mic_i(sizeof(mic_f)));
      const mic_i nearZoffset16 = select(rdir16.z >= 0.0f,mic_i(0),mic_i(sizeof(mic_f)));
#endif

      long rayIndex = -1;
      while((rayIndex = bitscan64(rayIndex,toInt(m_valid))) != BITSCAN_NO_BIT_SET_64)	    
        {
	  stack_node[1] = root;
	  size_t sindex = 2;

	  const mic_f org_xyz  = loadAOS4to16f(rayIndex,ray16.org.x,ray16.org.y,ray16.org.z);
	  const mic_f dir_xyz  = loadAOS4to16f(rayIndex,ray16.dir.x,ray16.dir.y,ray16.dir.z);

	  const mic3f rdir(rdir16.x[rayIndex],rdir16.y[rayIndex],rdir16.z[rayIndex]);
	  const mic_f rdir_xyz = rcp(dir_xyz);
	  const mic_f org_rdir_xyz = rdir_xyz * org_xyz;

	  const mic3f org(swAAAA(org_xyz),swBBBB(org_xyz),swCCCC(org_xyz));
	  const mic3f org_rdir = rdir * org;	  
	  const mic_f tnear = broadcast1to16f(&ray16.tnear[rayIndex]);
	  mic_f       tfar  = broadcast1to16f(&ray16.tfar[rayIndex]);

	  const unsigned int leaf_mask = TMP_BVH16I_LEAF_MASK;

#if defined(NEAR_FAR_OPT_OCCLUDED)            
	      const unsigned int nearXoffset0 = nearXoffset16[rayIndex];
	      const unsigned int nearYoffset0 = nearYoffset16[rayIndex];
	      const unsigned int nearZoffset0 = nearZoffset16[rayIndex];

	      const unsigned int nearXoffset1 = nearXoffset0 ^ sizeof(mic_f);
	      const unsigned int nearYoffset1 = nearYoffset0 ^ sizeof(mic_f);
	      const unsigned int nearZoffset1 = nearZoffset0 ^ sizeof(mic_f);
#endif

	  while (1)
	    {
	      NodeRef curNode = stack_node[sindex-1];
	      sindex--;
            
	      while (1) 
		{
		  /* test if this is a leaf node */
		  if (unlikely(curNode.isLeaf(leaf_mask))) break;
		  STAT3(shadow.trav_nodes,1,1,1);

		  const BVH16i::Node* __restrict__ const bptr = (BVH16i::Node*)((char*)bvh16 + curNode.id());

		  prefetch<PFHINT_L1>(&bptr->min_x);
		  prefetch<PFHINT_L1>(&bptr->max_x);
		  prefetch<PFHINT_L1>(&bptr->min_y);
		  prefetch<PFHINT_L1>(&bptr->max_y);
		  prefetch<PFHINT_L1>(&bptr->min_z);
		  prefetch<PFHINT_L1>(&bptr->max_z);
		  prefetch<PFHINT_L1>(&bptr->child);

#if !defined(NEAR_FAR_OPT_OCCLUDED)
        
		  const mic_f min_x = bptr->min_x * rdir.x - org_rdir.x;
		  const mic_f max_x = bptr->max_x * rdir.x - org_rdir.x;

		  const mic_f min_y = bptr->min_y * rdir.y - org_rdir.y;
		  const mic_f max_y = bptr->max_y * rdir.y - org_rdir.y;

		  const mic_f min_z = bptr->min_z * rdir.z - org_rdir.z;
		  const mic_f max_z = bptr->max_z * rdir.z - org_rdir.z;

		  const mic_f nearX = min(min_x,max_x);
		  const mic_f farX  = max(min_x,max_x);

		  const mic_f nearY = min(min_y,max_y);
		  const mic_f farY  = max(min_y,max_y);

		  const mic_f nearZ = min(min_z,max_z);
		  const mic_f farZ  = max(min_z,max_z);

#else

		  const mic_f nearX = load16f((float*)((const char*)&bptr->min_x + (size_t)nearXoffset0)) * rdir.x - org_rdir.x;
		  const mic_f farX  = load16f((float*)((const char*)&bptr->min_x + (size_t)nearXoffset1)) * rdir.x - org_rdir.x;

		  const mic_f nearY = load16f((float*)((const char*)&bptr->min_y + (size_t)nearYoffset0)) * rdir.y - org_rdir.y;
		  const mic_f farY  = load16f((float*)((const char*)&bptr->min_y + (size_t)nearYoffset1)) * rdir.y - org_rdir.y;

		  const mic_f nearZ = load16f((float*)((const char*)&bptr->min_z + (size_t)nearZoffset0)) * rdir.z - org_rdir.z;
		  const mic_f farZ  = load16f((float*)((const char*)&bptr->min_z + (size_t)nearZoffset1)) * rdir.z - org_rdir.z;
		  
#endif

		  const mic_f near16 = max(max(nearX,nearY),max(nearZ,tnear));
		  const mic_f far16  = min(min(farX,farY),min(farZ,tfar));

		  sindex--;

		  curNode = stack_node[sindex]; // early pop of next node

		  const mic_m hitm = le(near16,far16);
		  const mic_f tNear_pos = select(hitm,near16,inf);


		  // const BVH16i::Node* __restrict__ const next = (BVH16i::Node*)((char*)bvh16 + curNode.id());

		  // prefetch<PFHINT_L2>(&next->min_x);
		  // prefetch<PFHINT_L2>(&next->max_x);
		  // prefetch<PFHINT_L2>(&next->min_y);
		  // prefetch<PFHINT_L2>(&next->max_y);
		  // prefetch<PFHINT_L2>(&next->min_z);
		  // prefetch<PFHINT_L2>(&next->max_z);

		  STAT3(shadow.trav_hit_boxes[countbits(hitm)],1,1,1);


		  /* if no child is hit, continue with early popped child */
		  if (unlikely(none(hitm))) continue;


		  sindex++;
        
		  const unsigned long hiti = toInt(hitm);
		  const unsigned long pos_first = bitscan64(hiti);
		  const unsigned long num_hitm = countbits(hiti); 
        
		  /* if a single child is hit, continue with that child */
		  curNode = bptr->child[pos_first];
		  if (likely(num_hitm == 1)) continue;
        
		  /* if two children are hit, push in correct order */
		  const unsigned long pos_second = bitscan64(pos_first,hiti);
		  if (likely(num_hitm == 2))
		    {
		      const unsigned int dist_first  = ((unsigned int*)&near16)[pos_first];
		      const unsigned int dist_second = ((unsigned int*)&near16)[pos_second];
		      const unsigned int node_first  = curNode;
		      const unsigned int node_second = bptr->child[pos_second];
          
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
#if 1
		    {
		      const unsigned int old_sindex = sindex;
		      const mic_i children = bptr->child;
		      compactustore16i(hitm,&stack_node[old_sindex],children);
		      sindex += num_hitm-1;
		      curNode = stack_node[sindex];
		      continue;
		    }
#else

		  /* continue with closest child and push all others */
		  const mic_f min_dist = set_min16(tNear_pos);

		  const unsigned int old_sindex = sindex;
		  sindex += countbits(hiti) - 1;
		  assert(sindex < 3*BVH4i::maxDepth+1);
		  const mic_i children = bptr->child;
        
		  const mic_m closest_child = eq(hitm,min_dist,near16);
		  const unsigned long closest_child_pos = bitscan64(closest_child);
		  const mic_m m_pos = andn(hitm,andn(closest_child,(mic_m)((unsigned int)closest_child - 1)));
		  curNode = bptr->child[closest_child_pos];
		  compactustore16i(m_pos,&stack_node[old_sindex],children);        
#endif
		}
	  
	    
#if defined(ENABLE_EXTENDED_LEAVES)
	      if (likely((curNode.isLeaf(BVH16I_EXTENDED_LEAF_MASK) == BVH16I_EXTENDED_LEAF_MASK)))
		{
		  curNode.id() ^= BVH16I_EXTENDED_LEAF_MASK;

		  const BVH4i::Node* __restrict__ const node4 = (BVH4i::Node*)((char*)bvh16 + curNode.id());

		  const float* __restrict const plower = (float*)node4->lower;
		  const float* __restrict const pupper = (float*)node4->upper;

		  prefetch<PFHINT_L1>((char*)node4 + 0);
		  prefetch<PFHINT_L1>((char*)node4 + 64);
        
		  /* intersect single ray with 4 bounding boxes */
		  const mic_f tLowerXYZ = load16f(plower) * rdir_xyz - org_rdir_xyz;
		  const mic_f tUpperXYZ = load16f(pupper) * rdir_xyz - org_rdir_xyz;
		  const mic_f tLower = mask_min(0x7777,tnear,tLowerXYZ,tUpperXYZ);
		  const mic_f tUpper = mask_max(0x7777,tfar,tLowerXYZ,tUpperXYZ);

		  const mic_f tNear = vreduce_max4(tLower);
		  const mic_f tFar  = vreduce_min4(tUpper);  
		  const mic_m hitm = le(0x8888,tNear,tFar);
		  const mic_f tNear_pos = select(hitm,tNear,inf);

		  STAT3(normal.trav_hit_boxes[countbits(hitm)],1,1,1);

		  /* if no child is hit, continue with early popped child */
		  if (unlikely(none(hitm))) continue;
        
		  const unsigned long hiti = toInt(hitm);
		  const unsigned long pos_first = bitscan64(hiti);
		  const unsigned long num_hitm = countbits(hiti); 
        
		  /* if a single child is hit, continue with that child */
		  curNode = ((unsigned int *)plower)[pos_first];
		  if (likely(num_hitm == 1)) 
		    {
		      
		    }
		  else if (likely(num_hitm == 2)) /* if two children are hit, push in correct order */
		    {
		      const unsigned long pos_second = bitscan64(pos_first,hiti);
		      const unsigned int dist_first  = ((unsigned int*)&tNear)[pos_first];
		      const unsigned int dist_second = ((unsigned int*)&tNear)[pos_second];
		      const unsigned int node_first  = curNode;
		      const unsigned int node_second = ((unsigned int*)plower)[pos_second];
          
		      if (dist_first <= dist_second)
			{
			  stack_node[sindex] = node_second;
			  sindex++;
			  assert(sindex < 3*BVH4i::maxDepth+1);
			}
		      else
			{
			  stack_node[sindex] = curNode;
			  curNode = node_second;
			  sindex++;
			  assert(sindex < 3*BVH4i::maxDepth+1);
			}
		    }
		  else
		    {        
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
		      compactustore16i(m_pos,&stack_node[old_sindex],plower_node);
		    }		  
		}
#endif

	      /* return if stack is empty */
	      if (unlikely(curNode == BVH16I_TERMINAL_TOKEN)) break;

	      STAT3(shadow.trav_leaves,1,1,1);
	      STAT3(shadow.trav_prims,4,4,4);

	      /* intersect one ray against four triangles */

	      //////////////////////////////////////////////////////////////////////////////////////////////////

	      const Triangle1* tptr  = (Triangle1*) curNode.leaf(accel);
	      prefetch<PFHINT_L1>(tptr + 3);
	      prefetch<PFHINT_L1>(tptr + 2);
	      prefetch<PFHINT_L1>(tptr + 1);
	      prefetch<PFHINT_L1>(tptr + 0); 

	      const mic_i and_mask = broadcast4to16i(zlc4);
	      
	      const mic_f v0 = gather_4f_zlc(and_mask,
					     (float*)&tptr[0].v0,
					     (float*)&tptr[1].v0,
					     (float*)&tptr[2].v0,
					     (float*)&tptr[3].v0);
	      
	      const mic_f v1 = gather_4f_zlc(and_mask,
					     (float*)&tptr[0].v1,
					     (float*)&tptr[1].v1,
					     (float*)&tptr[2].v1,
					     (float*)&tptr[3].v1);
	      
	      const mic_f v2 = gather_4f_zlc(and_mask,
					     (float*)&tptr[0].v2,
					     (float*)&tptr[1].v2,
					     (float*)&tptr[2].v2,
					     (float*)&tptr[3].v2);

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

	      const mic_m valid_u = ge((mic_m)m_init,u,zero);
	      const mic_m valid_v = ge(valid_u,v,zero);
	      const mic_m m_aperture = le(valid_v,u+v,mic_f::one()); 

	      const mic_f nom = ldot3_zxy(org,normal);
	      const mic_f t = rcp_den*nom;
	      if (unlikely(none(m_aperture))) continue;

	      const mic_m m_final  = lt(lt(m_aperture,tnear,t),t,tfar);

	      if (unlikely(any(m_final)))
		{
#if defined(__USE_RAY_MASK__)
		  const mic_i rayMask(ray16.mask[rayIndex]);
		  const mic_i triMask = swDDDD(gather16i_4i((int*)&tptr[0].v2,
							    (int*)&tptr[1].v2,
							    (int*)&tptr[2].v2,
							    (int*)&tptr[3].v2));
		  const mic_m m_ray_mask = (rayMask & triMask) != mic_i::zero();
		    
		 if ( any(m_final & m_ray_mask) )
#endif

		   {
		     STAT3(shadow.trav_prim_hits,1,1,1);
		     terminated |= toMask(mic_m::shift1[rayIndex]);
		     break;
		   }
		}
	      //////////////////////////////////////////////////////////////////////////////////////////////////

	    }


	  if (unlikely(all(toMask(terminated)))) break;
	}


      store16i(m_valid & toMask(terminated),&ray16.geomID,0);
    }

    void BVH16iIntersector1::intersect(BVH4i* bvh, Ray& ray) {}
    void BVH16iIntersector1::occluded (BVH4i* bvh, Ray& ray) {}
    
    DEFINE_INTERSECTOR16   (BVH16iTriangle1Intersector16SingleMoeller, BVH16iIntersector16Single);
    DEFINE_INTERSECTOR1    (BVH16iTriangle1Intersector1, BVH16iIntersector1);

  }
}
