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

#include "bvh4aos_triangle1_intersector1_ref.h"
#include "../geometry/triangle1.h"

#define STAT_COLLECTOR(x)

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
  
  void BVH4AOSTriangle1Intersector1Ref::intersect(const BVH4AOSTriangle1Intersector1Ref* This, Ray& ray) //const mic_m m_active, Ray &ray)  {
  {
    MIC_ALIGN unsigned int stack_node[QUAD_BVH_MAX_STACK_DEPTH];
    MIC_ALIGN float stack_near[QUAD_BVH_MAX_STACK_DEPTH];
    MIC_ALIGN float pack_near[16];
    
    const BVH4AOS* bvh = This->bvh;
    const Node      *const __restrict__ qbvh  = (const Node*) bvh->nodePtr(); //Renderer::scene->qbvh;
    const Triangle1 *const __restrict__ accel = (const Triangle1*) bvh->triPtr(); //Renderer::scene->accel;
    
    assert((unsigned long)stack_node % 64 == 0);
    assert((unsigned long)stack_near % 64 == 0);
    
    prefetch<PFHINT_L1EX>(&stack_node[0]);      
    prefetch<PFHINT_L1EX>(&stack_node[16]);      
    prefetch<PFHINT_L1EX>(&stack_near[0]);      
    prefetch<PFHINT_L1EX>(&stack_near[16]);  

    //const mic3f &origin    = ray.org;
    //const mic3f &direction = ray.dir;
    //const mic_f &min_dist  = ray.t_min;
    //const mic_f  max_dist  = ray.tfar;
    
    //const mic_f direction_x = sel(eqz(direction.x),mic_f::ulp(),direction.x);
    //const mic_f direction_y = sel(eqz(direction.y),mic_f::ulp(),direction.y);
    //const mic_f direction_z = sel(eqz(direction.z),mic_f::ulp(),direction.z);
    
    //const mic_f rdirx = rcp(direction_x);
    //const mic_f rdiry = rcp(direction_y);
    //const mic_f rdirz = rcp(direction_z);
  
    //ray.triID    = mic_i::minus_one();
    //ray.shaderID = mic_i::minus_one();
        
    //long rayIndex = -1;
    //const unsigned int m_index = toInt(m_active);
    stack_node[0] = QBVH_TERMINAL_TOKEN;
    
    //while((rayIndex = bsf64(rayIndex,m_index)) != MIC_NO_BIT_SET_64) {
    stack_node[1] = bvh->root; //qbvh[0].lower[0].d;
    unsigned int sindex = 2;
      
    const mic_f inf = mic_f::inf();
    const Vector3f ray_rdir     = rcp_safe(ray.dir);
    const mic_f org_aos      = loadAOS(ray.org.x,ray.org.y,ray.org.z); //SOAtoAOS_4f(rayIndex,origin_x,origin_y,origin_z);
    const mic_f dir_aos      = loadAOS(ray.dir.x,ray.dir.y,ray.dir.z); //SOAtoAOS_4f(rayIndex,direction_x,direction_y,direction_z);
    const mic_f rdir_aos     = loadAOS(ray_rdir.x,ray_rdir.y,ray_rdir.z); //SOAtoAOS_4f(rayIndex,rdirx,rdiry,rdirz);
    const mic_f org_rdir_aos = org_aos * rdir_aos;
    const mic_f min_dist_aos = upconv1f(&ray.tnear); //upconv1f(&min_dist[rayIndex]);
      
    mic_i triangleID = mic_i::minus_one();
    mic_f max_dist_aos = upconv1f(&ray.tfar); //upconv1f(&max_dist[rayIndex]);
      
    *(mic_f*)stack_near = inf;
      
    STAT_COLLECTOR(StatCollector::single.numChunks++; StatCollector::single.numRays++);
      
    while (1)
    {
      unsigned int curNode = stack_node[sindex-1];
      sindex--;
      
        STAT_COLLECTOR(StatCollector::single.numStackPops++);
        
        while (1) 
        {
          if (unlikely(qbvhLeaf(curNode,/*TMP_*/QBVH_LEAF_MASK) != 0)) break;
          
          mic_f near4_min = inf;
          STAT_COLLECTOR(StatCollector::single.numTraversalSteps++);
          STAT_COLLECTOR(StatCollector::single.numBoxIntersections++);
          
          const Node* __restrict__ const qptr = qbvhChildPtrNoMask(qbvh,curNode);
          //const Node* __restrict__ const qptr = qbvhChildPtr(qbvh,curNode);
                    
          const float* __restrict const b_min = (float*)qptr->lower;
          const float* __restrict const b_max = (float*)qptr->upper;
          const mic_f clipMinXYZ = upconv16f(b_min) * rdir_aos - org_rdir_aos;
          const mic_f clipMaxXYZ = upconv16f(b_max) * rdir_aos - org_rdir_aos;
          
          curNode = stack_node[sindex-1];

          mic_f clipMinXYZ_min = min_dist_aos;
          clipMinXYZ_min = mask_min(0x7777,clipMinXYZ_min,clipMinXYZ,clipMaxXYZ);
          mic_f clipMaxXYZ_max = max_dist_aos;
          clipMaxXYZ_max = mask_max(0x7777,clipMaxXYZ_max,clipMinXYZ,clipMaxXYZ);
          
          sindex--;
          
          const mic_f near4 = set_max4(clipMinXYZ_min);
          const mic_f far4  = set_min4(clipMaxXYZ_max);  
          
          
          const mic_m m_lr_hit = le(0x8888,near4,far4);
          
          
          store16f(pack_near,near4);
          
          STAT_COLLECTOR(if ((unsigned int)m_lr_hit & 0x8) StatCollector::single.histo[0]++);
          STAT_COLLECTOR(if ((unsigned int)m_lr_hit & 0x80) StatCollector::single.histo[1]++);
          STAT_COLLECTOR(if ((unsigned int)m_lr_hit & 0x800) StatCollector::single.histo[2]++);
          STAT_COLLECTOR(if ((unsigned int)m_lr_hit & 0x8000) StatCollector::single.histo[3]++);
          
          near4_min = sel(m_lr_hit,near4,near4_min);
          
          if (unlikely(ez_mask(m_lr_hit) != 0)) continue;
          
          const unsigned long i_lr_hit = toInt(m_lr_hit);
          const unsigned long pos_first = bsf64(i_lr_hit);
          const unsigned long num_m_lr_hit = countbits(i_lr_hit); 
	  
          sindex++;
          assert(sindex < QUAD_BVH_MAX_STACK_DEPTH);
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
              stack_node[sindex] = node_sec;
              ((unsigned int*)stack_near)[sindex] = dist_sec;                      
              sindex++;
              assert(sindex < QUAD_BVH_MAX_STACK_DEPTH);
              continue;
            }
            else
            {
              assert(pack_near[pos_sec] <= pack_near[pos_first]);
              stack_node[sindex] = curNode;
              ((unsigned int*)stack_near)[sindex] = dist_first;
              curNode = node_sec;
              sindex++;
              assert(sindex < QUAD_BVH_MAX_STACK_DEPTH);
              continue;
            }
          }
          
          const mic_f child_min_dist = set_min_lanes(near4_min);
          
          const unsigned int *__restrict const b_index = (unsigned int*)b_min;
          const unsigned int old_sindex = sindex;
          sindex += countbits(i_lr_hit) - 1;
          assert(sindex < QUAD_BVH_MAX_STACK_DEPTH);
          
          const mic_m m_child_min_dist = eq(m_lr_hit,child_min_dist,near4);
          const unsigned long pos = bsf64(m_child_min_dist);
          //const mic_m m_pos = m_andnot(m_lr_hit,(mic_m)((unsigned int)1 << pos));
          
          
          
          
          const mic_m m_pos = andn(m_lr_hit,andn(m_child_min_dist,(mic_m)((unsigned int)m_child_min_dist - 1)));
          const mic_i b_min_node = load16i((int*)b_min);
          
          curNode = ((unsigned int*)b_min)[pos];
          
          compactustore16i(m_pos,&stack_node[old_sindex],b_min_node);
          compactustore16f(m_pos,&stack_near[old_sindex],near4);
          
        }
        
        if (unlikely(curNode == QBVH_TERMINAL_TOKEN)) { break; } // === stack empty ===
        
        STAT_COLLECTOR(StatCollector::single.numLeafIntersections++);
        
        unsigned int itemOffset    = qbvhItemOffset(curNode);
        const Triangle1  *__restrict__ tptr    = (Triangle1*)((char*)accel + itemOffset);
        
        prefetch<PFHINT_L1>(tptr + 3);
        prefetch<PFHINT_L1>(tptr + 2);
        prefetch<PFHINT_L1>(tptr + 1);
        prefetch<PFHINT_L1>(tptr + 0);
        
        
        const unsigned int items = qbvhItems(curNode);
        
        assert(items > 0);
        assert(items <= 4);
        
        unsigned int triID = qbvhItemOffsetToID(itemOffset);
        
        const mic_f zero = mic_f::zero();
        
        mic_i triID4 = upconv1i((const int*)&triID);
        
        mic_f u4 = undefined(); 
        mic_f v4 = undefined(); 
        
        mic_i triangleID4 = cast_to_mic_i(undefined());
        mic_m m_update = 0;
        const mic_i and_mask = mic_i::zlc4();
        prefetch<PFHINT_NT>(&triangleID4);      
        
        STAT_COLLECTOR(StatCollector::single.numAccessedPrimitives+=4);
        STAT_COLLECTOR(StatCollector::single.numPrimitiveIntersections++);
        
        
        /*const mic_f v0 = gather_4f_zlc(and_mask,
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
                                       (float*)&tptr[3].v2);*/

        const mic_f v0 = gather_4f_zlc(and_mask,(float*)&tptr[0].v0,(float*)&tptr[1].v0,(float*)&tptr[2].v0,(float*)&tptr[3].v0);
        const mic_f e1 = gather_4f_zlc(and_mask,(float*)&tptr[0].e1,(float*)&tptr[1].e1,(float*)&tptr[2].e1,(float*)&tptr[3].e1);
        const mic_f e2 = gather_4f_zlc(and_mask,(float*)&tptr[0].e2,(float*)&tptr[1].e2,(float*)&tptr[2].e2,(float*)&tptr[3].e2);
        //const mic_f Ng = gather_4f_zlc(and_mask,(float*)&tptr[0].Ng,(float*)&tptr[1].Ng,(float*)&tptr[2].Ng,(float*)&tptr[3].Ng);
        const mic_f normal = lcross_zxy(e2,e1);
                
        //const mic_f e1 = v1 - v0;
        //const mic_f e2 = v2 - v0;	     
        //const mic_f normal = lcross_zxy(e1,e2);

        const mic_f org = v0 - org_aos;
        const mic_f nom = ldot3_zxy(org,normal);
        const mic_m m_sign = lz(nom);
        
#if defined(BACK_FACE_CULLING)    
        if (unlikely(ez(m_sign) != 0)) continue;    
#endif
        
        const mic_f odzxy = msubr231(dir_aos * swizzle(org,_MM_SWIZ_REG_DACB),org,swizzle(dir_aos,_MM_SWIZ_REG_DACB));	      
        
        mic_m den_gt_0;
        const mic_f den = lz_ldot3_zxy(dir_aos,normal,den_gt_0);	      
        den_gt_0 = ~den_gt_0;
        
        const mic_f u = ldot3_zxy(e2,odzxy); 
        const mic_f v = ldot3_zxy(e1,odzxy); 
        
        
#if !defined(BACK_FACE_CULLING)    
        
        const mic_m valid_u = (den_gt_0^ge(u,zero));
        const mic_m valid_v = (den_gt_0^ge(v,zero));
        const mic_m valid_uv = den_gt_0^ge(add_neg(u,v),den); 
        const mic_m m_aperture = valid_u & valid_v & valid_uv;
        
#else
        const mic_m den_gt_0 = lt(m_sign,den,zero);
        const mic_m valid_u = ge(den_gt_0,u,zero); 
        const mic_m valid_v = ge(valid_u,v,zero); 
        const mic_m m_aperture = ge(valid_v,add_neg(u,v),den); 
#endif
        
        const mic_f rcp_den = rcp(den);
        
        if ( unlikely(ez_mask(m_aperture) != 0) ) continue;
        
        
        const mic_f t = rcp_den*nom;
	
        
        const mic_m m_final  = lt(lt(m_aperture ,min_dist_aos,t),t,max_dist_aos);
        const mic_i i_triID = triID4 + mic_i::addtriID4();
        
        max_dist_aos  = sel(m_final,t,max_dist_aos);
        
        triangleID4 = sel(m_final,i_triID,triangleID4); 
        //const mic_f rden = sel(den_gt_0,rcp_den,-rcp_den);
        
        const mic_f du = abs(u*rcp_den);		      
        const mic_f dv = abs(v*rcp_den);
        u4 = sel(m_final,du,u4); 
        v4 = sel(m_final,dv,v4);
        
        prefetch<PFHINT_L1EX>(&ray);
        //prefetch<PFHINT_L1EX>(&ray.tfar);   //prefetch<PFHINT_L1EX>(&ray.tfar);      
        //prefetch<PFHINT_L1EX>(&ray.id0); //prefetch<PFHINT_L1EX>(&ray.triID);      
        //prefetch<PFHINT_L1EX>(&ray.id1); //prefetch<PFHINT_L1EX>(&ray.shaderID);      
        
        m_update |= m_final;
        // ================ Update hit and compact stack ============
        if (unlikely(nz_mask(m_update) != 0))
        {
          //STAT_COLLECTOR(StatCollector::single.debug1++);
          
          //prefetch<PFHINT_L1EX>(&ray.u);      
          //prefetch<PFHINT_L1EX>(&ray.v);      
          //prefetch<PFHINT_L1EX>(&ray.gNormal[0]);      
          //prefetch<PFHINT_L1EX>(&ray.gNormal[1]);      
          //prefetch<PFHINT_L1EX>(&ray.gNormal[2]);      
          
          
          const mic_f min_dist = set_min16(max_dist_aos);
          
          const mic_m m_dist = eq(min_dist,max_dist_aos);
          
          max_dist_aos = min_dist;
          
          const unsigned long vecIndex = bsf64(m_dist);
          const unsigned long triIndex = triangleID4[vecIndex];
          const Triangle1  *__restrict__ tri_ptr = &accel[triIndex];
          
          const mic_m m_final = m_dist^(m_dist & (mic_m)((unsigned long)m_dist - 1));
          
          assert(countbits(m_final) == 1);
	  
          const mic_f gnormalz = mic_f(normal[vecIndex+0]);
          const mic_f gnormalx = mic_f(normal[vecIndex+1]);
          const mic_f gnormaly = mic_f(normal[vecIndex+2]);
          
          compactustore16f_low(m_final,&ray.tfar,min_dist); //compactustore16f_low(m_final,&ray.tfar[rayIndex],min_dist);
          //compactustore16i_low(m_final,&ray.triID[rayIndex],triangleID4);
          
          compactustore16f_low(m_final,&ray.u,u4); //compactustore16f_low(m_final,&ray.u[rayIndex],u4);
          compactustore16f_low(m_final,&ray.v,v4); //compactustore16f_low(m_final,&ray.v[rayIndex],v4);
          
          compactustore16f_low(m_final,&ray.Ng.z,gnormalz);
          compactustore16f_low(m_final,&ray.Ng.x,gnormalx);
          compactustore16f_low(m_final,&ray.Ng.y,gnormaly);
          
          //ray.shaderID[(unsigned long)rayIndex] = tri_ptr->shaderID;
          ray.id0 = tri_ptr->e1.a;
          ray.id1 = tri_ptr->e2.a;
          
          // === only after max_dist_aos is updated, do a stack compaction ===
          
          if (likely(sindex >= 2))
          {
            //STAT_COLLECTOR(StatCollector::single.debug0++);
            
            unsigned int old_sindex = sindex;
            if (likely(sindex < 16))
            {
              //const unsigned int m_num_stack = (((unsigned int)1 << sindex) - 1);
              const unsigned int m_num_stack = mic_i::shift1(sindex) - 1;
              
              const mic_m m_num_stack_low  = toMask(m_num_stack);
              const mic_f snear_low  = upconv16f(stack_near + 0);
              const mic_i snode_low  = upconv16i((int*)stack_node + 0);
              const mic_m m_stack_compact_low  = le(m_num_stack_low,snear_low,max_dist_aos) | (mic_m)1;
              
              compactustore16f_low(m_stack_compact_low,stack_near + 0,snear_low);
              compactustore16i_low(m_stack_compact_low,(int*)stack_node + 0,snode_low);
              sindex = countbits(m_stack_compact_low);
              assert(sindex < 32);
              
            }
            else
            {
              //const unsigned int m_num_stack = (((unsigned int)1 << sindex) - 1);
              const unsigned int m_num_stack = mic_i::shift1(sindex) - 1;
              
              const mic_m m_num_stack_low  = toMask(m_num_stack);
              const mic_m m_num_stack_high = toMask(m_num_stack >> 16);
              
              const mic_f snear_low  = upconv16f(stack_near + 0);
              const mic_f snear_high = upconv16f(stack_near + 16);
              
              const mic_i snode_low  = upconv16i((int*)stack_node + 0);
              const mic_i snode_high = upconv16i((int*)stack_node + 16);
              
              const mic_m m_stack_compact_low  = le(m_num_stack_low,snear_low,max_dist_aos) | (mic_m)1;
              const mic_m m_stack_compact_high = le(m_num_stack_high,snear_high,max_dist_aos);
              
              compactustore16f(m_stack_compact_low,stack_near + 0,snear_low);
              compactustore16i(m_stack_compact_low,(int*)stack_node + 0,snode_low);
              
              compactustore16f(m_stack_compact_high,stack_near + countbits(m_stack_compact_low),snear_high);
              compactustore16i(m_stack_compact_high,(int*)stack_node + countbits(m_stack_compact_low),snode_high);
              
              sindex = countbits(m_stack_compact_low) + countbits(m_stack_compact_high);
              assert(sindex < 32);
            }
          }
        }
      }
      //}
  }
  
  bool BVH4AOSTriangle1Intersector1Ref::occluded(const BVH4AOSTriangle1Intersector1Ref* This, Ray& ray) //const mic_m m_active, Ray &ray)
  {
    MIC_ALIGN unsigned int stack_node[QUAD_BVH_MAX_STACK_DEPTH];
    MIC_ALIGN float pack_near[16];

    const BVH4AOS* bvh = This->bvh;
    const Node      *const __restrict__ qbvh  = (const Node*) bvh->nodePtr(); //Renderer::scene->qbvh;
    const Triangle1 *const __restrict__ accel = (const Triangle1*) bvh->triPtr(); //Renderer::scene->accel;
    
    prefetch<PFHINT_L1EX>(&stack_node[0]);      
    prefetch<PFHINT_L1EX>(&stack_node[16]);      
    
    //const mic3f &origin    = ray.org;
    //const mic3f &direction = ray.dir;
    //const mic_f &min_dist  = ray.t_min;
    //const mic_f  max_dist  = ray.tfar;
    
    //const mic_f direction_x = sel(eqz(direction.x),mic_f::ulp(),direction.x);
    //const mic_f direction_y = sel(eqz(direction.y),mic_f::ulp(),direction.y);
    //const mic_f direction_z = sel(eqz(direction.z),mic_f::ulp(),direction.z);
    
    //const mic_f rdirx = rcp(direction_x);
    //const mic_f rdiry = rcp(direction_y);
    //const mic_f rdirz = rcp(direction_z);
    
    //mic_i not_occluded = mic_i::minus_one();
    
    //long rayIndex = -1;
    //const unsigned int m_index = toInt(m_active);
    stack_node[0] = QBVH_TERMINAL_TOKEN;
    
    //while((rayIndex = bsf64(rayIndex,m_index)) != MIC_NO_BIT_SET_64) 
    //{
    stack_node[1] = bvh->root; //qbvh[0].lower[0].d;
    unsigned int sindex = 2;
    
    const mic_f inf = mic_f::inf();
    const Vector3f ray_rdir     = rcp_safe(ray.dir);
    const mic_f org_aos      = loadAOS(ray.org.x,ray.org.y,ray.org.z); //SOAtoAOS_4f(rayIndex,origin_x,origin_y,origin_z);
    const mic_f dir_aos      = loadAOS(ray.dir.x,ray.dir.y,ray.dir.z); //SOAtoAOS_4f(rayIndex,direction_x,direction_y,direction_z);
    const mic_f rdir_aos     = loadAOS(ray_rdir.x,ray_rdir.y,ray_rdir.z); //SOAtoAOS_4f(rayIndex,rdirx,rdiry,rdirz);
    const mic_f org_rdir_aos = org_aos * rdir_aos;
    const mic_f min_dist_aos = upconv1f(&ray.tnear); //upconv1f(&min_dist[rayIndex]);
    const mic_f max_dist_aos = upconv1f(&ray.tfar); //upconv1f(&max_dist[rayIndex]);
    
    STAT_COLLECTOR(StatCollector::single.numChunks++; StatCollector::single.numRays++);
    
    const mic_m m_lastInLane = 0x8888;
    const unsigned int tmp_mask = /*TMP_*/QBVH_LEAF_MASK;
    while (1)
    {
      unsigned int curNode = stack_node[sindex-1];
      
      sindex--;
      
      STAT_COLLECTOR(StatCollector::single.numStackPops++);
      
      
      while (1) 
      {
        if (unlikely(qbvhLeaf(curNode,tmp_mask) != 0)) break;
        
        const Node* __restrict__ const qptr = qbvhChildPtrNoMask(qbvh,curNode);
        
        const float* __restrict const b_min = (float*)qptr->lower;
        const float* __restrict const b_max = (float*)qptr->upper;
        
        //prefetch<PFHINT_L1>(b_max); 
        //prefetch<PFHINT_L1>(b_min); 
        
        mic_f near4_min = inf;
        STAT_COLLECTOR(StatCollector::single.numTraversalSteps++);
        STAT_COLLECTOR(StatCollector::single.numBoxIntersections++);
        
        
        const mic_f clipMinXYZ = upconv16f(b_min) * rdir_aos - org_rdir_aos;
        const mic_f clipMaxXYZ = upconv16f(b_max) * rdir_aos - org_rdir_aos;
        
        curNode = stack_node[sindex-1];
        
        const mic_m m_threeLines = ~m_lastInLane;
        mic_f clipMinXYZ_min = min_dist_aos;
        clipMinXYZ_min = mask_min(m_threeLines,clipMinXYZ_min,clipMinXYZ,clipMaxXYZ);
        mic_f clipMaxXYZ_max = max_dist_aos;
        clipMaxXYZ_max = mask_max(m_threeLines,clipMaxXYZ_max,clipMinXYZ,clipMaxXYZ);
        
        sindex--;
        
        
        
        
        const mic_f near4 = set_max4(clipMinXYZ_min);
        const mic_f far4  = set_min4(clipMaxXYZ_max);  
        
        
        
        const mic_m m_lr_hit = le(m_lastInLane,near4,far4); 
        
        store16f(pack_near,near4);
        
        
        near4_min = sel(m_lr_hit,near4,near4_min);
        
        if (unlikely(ez_mask(m_lr_hit) != 0)) continue;
        const unsigned long i_lr_hit = toInt(m_lr_hit);
        
        const unsigned long pos_first = bsf64(i_lr_hit);
        const unsigned long num_m_lr_hit = countbits(i_lr_hit); 
	
        sindex++;
        assert(sindex < QUAD_BVH_MAX_STACK_DEPTH);
        curNode = ((unsigned int*)b_min)[pos_first];
        
        
        if (likely(num_m_lr_hit == 1)) continue;
#if 1
        
        
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
            stack_node[sindex] = node_sec;
            sindex++;
            assert(sindex < QUAD_BVH_MAX_STACK_DEPTH);
            
            continue;
          }
          else
          {
            assert(pack_near[pos_sec] <= pack_near[pos_first]);
            stack_node[sindex] = curNode;
            curNode = node_sec;
            sindex++;
            assert(sindex < QUAD_BVH_MAX_STACK_DEPTH);
            
            continue;
          }
        }
#endif
        
        const mic_f child_min_dist = set_min_lanes(near4_min);
        
        
        const unsigned int *__restrict const b_index = (unsigned int*)b_min;
        
        const unsigned int old_sindex = sindex;
        sindex += countbits(i_lr_hit) - 1;
        assert(sindex < QUAD_BVH_MAX_STACK_DEPTH);
        
        const mic_m m_child_min_dist = eq(m_lr_hit,child_min_dist,near4);
        const unsigned long pos = bsf64(m_child_min_dist);
        
        const mic_m m_pos = andn(m_lr_hit,andn(m_child_min_dist,(mic_m)((unsigned int)m_child_min_dist - 1)));
        const mic_i b_min_node = load16i((int*)b_min);
        
        curNode = ((unsigned int*)b_min)[pos];
        
        
        
        compactustore16i(m_pos,&stack_node[old_sindex],b_min_node);
      }
      
      if (unlikely(curNode == QBVH_TERMINAL_TOKEN)) { break; } // === stack empty ===
      
      STAT_COLLECTOR(StatCollector::single.numLeafIntersections++);
      
      unsigned int itemOffset    = qbvhItemOffset(curNode);
      const Triangle1  *__restrict__ tptr    = (Triangle1*)((char*)accel + itemOffset);
      
      prefetch<PFHINT_NT>(tptr + 3);
      prefetch<PFHINT_NT>(tptr + 2);
      prefetch<PFHINT_NT>(tptr + 1);
      prefetch<PFHINT_NT>(tptr + 0);
      
      
      const unsigned int items = qbvhItems(curNode);
      
      assert(items > 0);
      assert(items <= 4);
      
      const mic_f zero = mic_f::zero();
      
      //mic_m m_update = 0;
      const mic_i and_mask = mic_i::zlc4();
      
      {
        STAT_COLLECTOR(StatCollector::single.numAccessedPrimitives+=4);
        STAT_COLLECTOR(StatCollector::single.numPrimitiveIntersections++);
        
        
        /*const mic_f v0 = gather_4f_zlc(and_mask,
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
                                       (float*)&tptr[3].v2);*/

        const mic_f v0 = gather_4f_zlc(and_mask,(float*)&tptr[0].v0,(float*)&tptr[1].v0,(float*)&tptr[2].v0,(float*)&tptr[3].v0);
        const mic_f e1 = gather_4f_zlc(and_mask,(float*)&tptr[0].e1,(float*)&tptr[1].e1,(float*)&tptr[2].e1,(float*)&tptr[3].e1);
        const mic_f e2 = gather_4f_zlc(and_mask,(float*)&tptr[0].e2,(float*)&tptr[1].e2,(float*)&tptr[2].e2,(float*)&tptr[3].e2);
        //const mic_f Ng = gather_4f_zlc(and_mask,(float*)&tptr[0].Ng,(float*)&tptr[1].Ng,(float*)&tptr[2].Ng,(float*)&tptr[3].Ng);
        const mic_f normal = lcross_zxy(e2,e1);
                
        //const mic_f e1 = v1 - v0;
        //const mic_f e2 = v2 - v0;	     
        //const mic_f normal = lcross_zxy(e1,e2);

        const mic_f org = v0 - org_aos;
        const mic_f nom = ldot3_zxy(org,normal);
        const mic_m m_sign = lz(nom);
        
#if defined(BACK_FACE_CULLING)    
        if (unlikely(ez(m_sign) != 0)) continue;    
#endif
        
        const mic_f odzxy = msubr231(dir_aos * swizzle(org,_MM_SWIZ_REG_DACB),org,swizzle(dir_aos,_MM_SWIZ_REG_DACB));	      
        
        mic_m den_gt_0;
        const mic_f den = lz_ldot3_zxy(dir_aos,normal,den_gt_0);	      
        den_gt_0 = ~den_gt_0;
        
        const mic_f u = ldot3_zxy(e2,odzxy); 
        const mic_f v = ldot3_zxy(e1,odzxy); 
	
        
        
#if !defined(BACK_FACE_CULLING)    
        
        const mic_m valid_u = (den_gt_0^ge(u,zero));
        const mic_m valid_v = (den_gt_0^ge(v,zero));
        const mic_m valid_uv = den_gt_0^ge(add_neg(u,v),den); 
        const mic_m m_aperture = valid_u & valid_v & valid_uv;
        
#else
        const mic_m den_gt_0 = lt(m_sign,den,zero);
        const mic_m valid_u = ge(den_gt_0,u,zero); 
        const mic_m valid_v = ge(valid_u,v,zero); 
        const mic_m m_aperture = ge(valid_v,add_neg(u,v),den); 
#endif
        
        const mic_f rcp_den = rcp(den);
        const mic_f t = rcp_den*nom;	     
        
        if ( unlikely(ez_mask(m_aperture) != 0) ) continue;
        
        mic_m m_final  = lt(lt(m_aperture ,min_dist_aos,t),t,max_dist_aos);
        
#if defined(ENABLE_LIGHTSOURCES_MASK)
        
        const mic_i lsMask0 = mic_i(tptr[0].data[0]);
        const mic_i lsMask1 = mic_i(tptr[1].data[0]);
        const mic_i lsMask2 = mic_i(tptr[2].data[0]);
        const mic_i lsMask3 = mic_i(tptr[3].data[0]);
        mic_i lsMask = lsMask0;
        lsMask = sel(0xf0  ,lsMask1,lsMask);
        lsMask = sel(0xf00 ,lsMask2,lsMask);
        lsMask = sel(0xf000,lsMask3,lsMask);
        m_final = mask_test(m_final,lsMask,mic_i(ray.data[rayIndex]));
        
#endif	      
        if (unlikely(m_final != 0))
        {
          //not_occluded[rayIndex] = 0;
          return true;
          //break;
        }
        
      }
      
      
      
      
    }
    
    //}
    //return gez(m_active,not_occluded);
    return false;
  }

  void BVH4AOSTriangle1Intersector1RefRegister () {
    TriangleMesh::intersectors1.add("bvh4aos","triangle1","fast_ref","moeller",true ,BVH4AOSTriangle1Intersector1Ref::create);
  }
}

