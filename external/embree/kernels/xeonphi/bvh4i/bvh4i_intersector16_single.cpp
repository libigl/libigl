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

#include "bvh4i_intersector16_single.h"
#include "geometry/triangle1.h"
#include "geometry/filter.h"

namespace embree
{
  namespace isa
  {

    static unsigned int BVH4I_LEAF_MASK = BVH4i::leaf_mask; // needed due to compiler efficiency bug

    static __aligned(64) int zlc4[4] = {0xffffffff,0xffffffff,0xffffffff,0};

    void BVH4iIntersector16Single::intersect(mic_i* valid_i, BVH4i* bvh, Ray16& ray16)
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
            
	      while (1) 
		{
		  /* test if this is a leaf node */
		  if (unlikely(curNode.isLeaf(leaf_mask))) break;
		  STAT3(normal.trav_nodes,1,1,1);
        
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

		  STAT3(normal.trav_hit_boxes[countbits(hitm)],1,1,1);

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
	  
	    

	      /* return if stack is empty */
	      if (unlikely(curNode == BVH4i::invalidNode)) break;

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

	      mic_m m_final  = lt(lt(m_aperture,min_dist_xyz,t),t,max_dist_xyz);


#if defined(__USE_RAY_MASK__)
	      const mic_i rayMask(ray16.mask[rayIndex]);
	      const mic_i triMask = swDDDD(gather16i_4i_align(&tptr[0].v2,&tptr[1].v2,&tptr[2].v2,&tptr[3].v2));
	      const mic_m m_ray_mask = (rayMask & triMask) != mic_i::zero();
	      m_final &= m_ray_mask;	      
#endif

	      //////////////////////////////////////////////////////////////////////////////////////////////////

              /* did the ray hit one of the four triangles? */
	      if (unlikely(any(m_final)))
		{
		  /* intersection filter test */
#if defined(__INTERSECTION_FILTER__) 
		  mic_f org_max_dist_xyz = max_dist_xyz;

		  /* did the ray hit one of the four triangles? */
		  while (any(m_final)) 
		    {
		      max_dist_xyz  = select(m_final,t,org_max_dist_xyz);
		      const mic_f min_dist = vreduce_min(max_dist_xyz);
		      const mic_m m_dist = eq(min_dist,max_dist_xyz);
		      const size_t vecIndex = bitscan(toInt(m_dist));
		      const size_t triIndex = vecIndex >> 2;
		      const Triangle1  *__restrict__ tri_ptr = tptr + triIndex;
		      const mic_m m_tri = m_dist^(m_dist & (mic_m)((unsigned int)m_dist - 1));
		      const mic_f gnormalx = mic_f(tri_ptr->Ng.x);
		      const mic_f gnormaly = mic_f(tri_ptr->Ng.y);
		      const mic_f gnormalz = mic_f(tri_ptr->Ng.z);
		      const int geomID = tri_ptr->geomID();
		      const int primID = tri_ptr->primID();
                
		      Geometry* geom = ((Scene*)bvh->geometry)->get(geomID);
		      if (likely(!geom->hasIntersectionFilter16())) 
			{

			  compactustore16f_low(m_tri,&ray16.tfar[rayIndex],min_dist);
			  compactustore16f_low(m_tri,&ray16.u[rayIndex],u); 
			  compactustore16f_low(m_tri,&ray16.v[rayIndex],v); 
			  compactustore16f_low(m_tri,&ray16.Ng.x[rayIndex],gnormalx); 
			  compactustore16f_low(m_tri,&ray16.Ng.y[rayIndex],gnormaly); 
			  compactustore16f_low(m_tri,&ray16.Ng.z[rayIndex],gnormalz); 
			  ray16.geomID[rayIndex] = geomID;
			  ray16.primID[rayIndex] = primID;
			  max_dist_xyz = min_dist;
			  break;
			}
                
		      if (runIntersectionFilter16(geom,ray16,rayIndex,u,v,min_dist,gnormalx,gnormaly,gnormalz,m_tri,geomID,primID)) {
			max_dist_xyz = min_dist;
			break;
		      }
		      m_final ^= m_tri;
		    }
		  max_dist_xyz = ray16.tfar[rayIndex];
#else
		  STAT3(normal.trav_prim_hits,1,1,1);
                  max_dist_xyz  = select(m_final,t,max_dist_xyz);
		  const mic_f min_dist = vreduce_min(max_dist_xyz);
		  const mic_m m_dist = eq(min_dist,max_dist_xyz);
                  const size_t vecIndex = bitscan(toInt(m_dist));
		  const size_t triIndex = vecIndex >> 2;

		  const Triangle1  *__restrict__ tri_ptr = tptr + triIndex;

		  const mic_m m_tri = m_dist^(m_dist & (mic_m)((unsigned int)m_dist - 1));

		  const mic_f gnormalx = mic_f(tri_ptr->Ng.x);
		  const mic_f gnormaly = mic_f(tri_ptr->Ng.y);
		  const mic_f gnormalz = mic_f(tri_ptr->Ng.z);
                
		  prefetch<PFHINT_L1EX>(&ray16.tfar);  
		  prefetch<PFHINT_L1EX>(&ray16.u);
		  prefetch<PFHINT_L1EX>(&ray16.v);
		  prefetch<PFHINT_L1EX>(&ray16.Ng.x); 
		  prefetch<PFHINT_L1EX>(&ray16.Ng.y); 
		  prefetch<PFHINT_L1EX>(&ray16.Ng.z); 
		  prefetch<PFHINT_L1EX>(&ray16.geomID);
		  prefetch<PFHINT_L1EX>(&ray16.primID);

		  max_dist_xyz = min_dist;
		  
		  compactustore16f_low(m_tri,&ray16.tfar[rayIndex],min_dist);
		  compactustore16f_low(m_tri,&ray16.u[rayIndex],u); 
		  compactustore16f_low(m_tri,&ray16.v[rayIndex],v); 
		  compactustore16f_low(m_tri,&ray16.Ng.x[rayIndex],gnormalx); 
		  compactustore16f_low(m_tri,&ray16.Ng.y[rayIndex],gnormaly); 
		  compactustore16f_low(m_tri,&ray16.Ng.z[rayIndex],gnormalz); 

		  ray16.geomID[rayIndex] = tri_ptr->geomID();
		  ray16.primID[rayIndex] = tri_ptr->primID();
#endif
		  /* compact the stack if size if stack >= 2 */
		  // ------------------------
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
	      // ------------------------
	    }	  
	}
    }
    
    void BVH4iIntersector16Single::occluded(mic_i* valid_i, BVH4i* bvh, Ray16& ray16)
    {
      /* near and node stack */
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

      /* setup */
      const mic_m m_valid     = *(mic_i*)valid_i != mic_i(0);
      const mic3f rdir16      = rcp_safe(ray16.dir);
      unsigned int terminated = toInt(!m_valid);
      const mic_f inf         = mic_f(pos_inf);
      const mic_f zero        = mic_f::zero();

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

	  const unsigned int leaf_mask = BVH4I_LEAF_MASK;

	  while (1)
	    {
	      NodeRef curNode = stack_node[sindex-1];
	      sindex--;
            
	      while (1) 
		{
		  /* test if this is a leaf node */
		  if (unlikely(curNode.isLeaf(leaf_mask))) break;
		  STAT3(shadow.trav_nodes,1,1,1);

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

		  const Node* __restrict__ const next = curNode.node(nodes);
		  prefetch<PFHINT_L2>((char*)next + 0);
		  prefetch<PFHINT_L2>((char*)next + 64);

		  sindex--;
		  const mic_f tNear = vreduce_max4(tLower);
		  const mic_f tFar  = vreduce_min4(tUpper);  
		  const mic_m hitm = le(0x8888,tNear,tFar);
		  const mic_f tNear_pos = select(hitm,tNear,inf);

		  curNode = stack_node[sindex]; // early pop of next node

		  STAT3(shadow.trav_hit_boxes[countbits(hitm)],1,1,1);

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
	  
	    

	      /* return if stack is empty */
	      if (unlikely(curNode == BVH4i::invalidNode)) break;

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

	      mic_m m_final  = lt(lt(m_aperture,min_dist_xyz,t),t,max_dist_xyz);

#if defined(__USE_RAY_MASK__)
	      const mic_i rayMask(ray16.mask[rayIndex]);
	      const mic_i triMask = swDDDD(gather16i_4i_align(&tptr[0].v2,&tptr[1].v2,&tptr[2].v2,&tptr[3].v2));
	      const mic_m m_ray_mask = (rayMask & triMask) != mic_i::zero();
	      m_final &= m_ray_mask;	      
#endif

#if defined(__INTERSECTION_FILTER__) 
              
              /* did the ray hit one of the four triangles? */
              while (any(m_final)) 
		{
		  const mic_f temp_t  = select(m_final,t,max_dist_xyz);
		  const mic_f min_dist = vreduce_min(temp_t);
		  const mic_m m_dist = eq(min_dist,temp_t);
		  const size_t vecIndex = bitscan(toInt(m_dist));
		  const size_t triIndex = vecIndex >> 2;
		  const Triangle1  *__restrict__ tri_ptr = tptr + triIndex;
		  const mic_m m_tri = m_dist^(m_dist & (mic_m)((unsigned int)m_dist - 1));
		  const mic_f gnormalx = mic_f(tri_ptr->Ng.x);
		  const mic_f gnormaly = mic_f(tri_ptr->Ng.y);
		  const mic_f gnormalz = mic_f(tri_ptr->Ng.z);
		  const int geomID = tri_ptr->geomID();
		  const int primID = tri_ptr->primID();                
		  Geometry* geom = ((Scene*)bvh->geometry)->get(geomID);
		  if (likely(!geom->hasOcclusionFilter16())) break;
                
		  if (runOcclusionFilter16(geom,ray16,rayIndex,u,v,min_dist,gnormalx,gnormaly,gnormalz,m_tri,geomID,primID)) 
		    break;

		  m_final ^= m_tri; /* clear bit */
		}
#endif
	      if (unlikely(any(m_final)))
		{
		  STAT3(shadow.trav_prim_hits,1,1,1);
		  terminated |= mic_m::shift1[rayIndex];
		  break;
		}
	      //////////////////////////////////////////////////////////////////////////////////////////////////

	    }


	  if (unlikely(all(toMask(terminated)))) break;
	}


      store16i(m_valid & toMask(terminated),&ray16.geomID,0);
    }
    
    DEFINE_INTERSECTOR16    (BVH4iTriangle1Intersector16SingleMoeller, BVH4iIntersector16Single);

  }
}
