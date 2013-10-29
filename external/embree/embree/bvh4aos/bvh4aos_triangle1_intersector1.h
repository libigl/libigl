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

#ifndef __EMBREE_BVH4AOS_TRIANGLE1_INTERSECTOR1_H__
#define __EMBREE_BVH4AOS_TRIANGLE1_INTERSECTOR1_H__

#include "bvh4aos.h"
#include "../include/intersector1.h"
#include "../geometry/triangles.h"

namespace embree
{
  /*! BVH4AOS Traverser. Single ray traversal implementation for a Quad BVH. */
  class BVH4AOSTriangle1Intersector1 : public Intersector1
  {
    /* shortcuts for frequently used types */
    typedef typename BVH4AOS::NodeRef NodeRef;
    typedef typename BVH4AOS::Node Node;
    
  public:
    BVH4AOSTriangle1Intersector1 (const BVH4AOS* bvh) 
      : Intersector1((intersectFunc)intersect,(occludedFunc)occluded), bvh(bvh) {}

    static Intersector1* create(const Accel* bvh) { 
      return new BVH4AOSTriangle1Intersector1((const BVH4AOS*)bvh); 
    }

    static void intersect(const BVH4AOSTriangle1Intersector1* This, Ray& ray);
    static bool occluded (const BVH4AOSTriangle1Intersector1* This, Ray& ray);
   
  public:
    static __forceinline void intersect1(const BVH4AOS* bvh, 
                                  const Node*     const nodes,
                                  const Triangle1 *const triangles,
                                  unsigned rootNode, long idx, 
                                  const mic_f org_xyz, const mic_f dir_xyz, 
                                  const mic_f rdir_xyz, const mic_f org_rdir_xyz, 
                                  const mic_f min_dist_xyz, mic_f max_dist_xyz, 
                                  Ray16& ray)
  {
    /* allocate stack and push root node */
    MIC_ALIGN unsigned stack_node[3*BVH4AOS::maxDepth+1];
    MIC_ALIGN float    stack_near[3*BVH4AOS::maxDepth+1];
    prefetch<PFHINT_L1EX>(&stack_node[0]);      
    prefetch<PFHINT_L1EX>(&stack_node[16]);      
    prefetch<PFHINT_L1EX>(&stack_near[0]);      
    prefetch<PFHINT_L1EX>(&stack_near[16]);  
    stack_node[0] = BVH4AOS::invalid_node;
    stack_node[1] = rootNode;
    *(mic_f*)stack_near = mic_f::inf();
    unsigned sindex = 2;
    const unsigned leaf_mask = BVH4AOS_LEAF_MASK;

    while (1)
    {
      unsigned cur = stack_node[sindex-1];
      sindex--;
            
      while (1) 
      {
        /* test if this is a leaf node */
        if (unlikely(cur & leaf_mask)) 
          break;
        
        const Node* __restrict__ const node = (Node*)((char*)nodes + (unsigned long)cur);
        const float* __restrict const plower = (float*)node->lower;
        const float* __restrict const pupper = (float*)node->upper;
        prefetch<PFHINT_L1>((char*)nodes + (unsigned long)cur);
        prefetch<PFHINT_L1>((char*)nodes + (unsigned long)cur + 64);
        
        /* intersect single ray with 4 bounding boxes */
        const mic_f tLowerXYZ = upconv16f(plower) * rdir_xyz - org_rdir_xyz;
        const mic_f tUpperXYZ = upconv16f(pupper) * rdir_xyz - org_rdir_xyz;
        cur = stack_node[sindex-1]; // early pop of next node
        const mic_f tLower = mask_min(0x7777,min_dist_xyz,tLowerXYZ,tUpperXYZ);
        const mic_f tUpper = mask_max(0x7777,max_dist_xyz,tLowerXYZ,tUpperXYZ);
        prefetch<PFHINT_L2>((char*)nodes + (unsigned long)cur);
        prefetch<PFHINT_L2>((char*)nodes + (unsigned long)cur + 64);
        sindex--;
        const mic_f tNear = set_max4(tLower);
        const mic_f tFar  = set_min4(tUpper);  
        const mic_m hitm = le(0x8888,tNear,tFar);
        
        /* if no child is hit, continue with early popped child */
        if (unlikely(ez_mask(hitm) != 0)) continue;
        sindex++;
        
        const unsigned long hiti = toInt(hitm);
        const unsigned long pos_first = bsf64(hiti);
        const unsigned long num_hitm = countbits(hiti); 
        
        /* if a single child is hit, continue with that child */
        cur = ((unsigned*)plower)[pos_first];
        if (likely(num_hitm == 1)) continue;
        
        /* if two children are hit, push in correct order */
        const unsigned long pos_second = bsf64(pos_first,hiti);
        if (likely(num_hitm == 2))
        {
          const unsigned dist_first = ((unsigned*)&tNear)[pos_first];
          const unsigned dist_second   = ((unsigned*)&tNear)[pos_second];
          const unsigned node_first = cur;
          const unsigned node_second   = ((unsigned*)plower)[pos_second];
          
          if (dist_first <= dist_second)
          {
            stack_node[sindex] = node_second;
            ((unsigned*)stack_near)[sindex] = dist_second;                      
            sindex++;
            assert(sindex < 3*BVH4AOS::maxDepth+1);
            continue;
          }
          else
          {
            stack_node[sindex] = cur;
            ((unsigned*)stack_near)[sindex] = dist_first;
            cur = node_second;
            sindex++;
            assert(sindex < 3*BVH4AOS::maxDepth+1);
            continue;
          }
        }
        
        /* continue with closest child and push all others */
        const mic_f min_dist = set_min_lanes(sel(hitm,tNear,inf));
        const unsigned old_sindex = sindex;
        sindex += countbits(hiti) - 1;
        assert(sindex < 3*BVH4AOS::maxDepth+1);
        
        const mic_m closest_child = eq(hitm,min_dist,tNear);
        const unsigned long closest_child_pos = bsf64(closest_child);
        const mic_m m_pos = andn(hitm,andn(closest_child,(mic_m)((unsigned)closest_child - 1)));
        const mic_i plower_node = load16i((int*)plower);
        cur = ((unsigned*)plower)[closest_child_pos];
        compactustore16i(m_pos,&stack_node[old_sindex],plower_node);
        compactustore16f(m_pos,&stack_near[old_sindex],tNear);
      }
      
      /* return if stack is empty */
      if (unlikely(cur == BVH4AOS_INVALID_NODE)) 
        break;
      
      /* decode leaf node */
      const unsigned itemOffset = BVH4AOS::itemOffset(cur);
      const unsigned items      = BVH4AOS::itemCount (cur);
      const Triangle1* triangle = (Triangle1*)((char*)triangles + itemOffset);

      /* intersect triangles */
      prefetch<PFHINT_L1>(triangle + 3);
      prefetch<PFHINT_L1>(triangle + 2);
      prefetch<PFHINT_L1>(triangle + 1);
      prefetch<PFHINT_L1>(triangle + 0);
      assert(items > 0);
      assert(items <= 4);
      
      /* calculate edges, geometry normal, and determinant */
      const mic_i and_mask = mic_i::zlc4();
      const mic_f p0 = gather_4f_zlc(and_mask,(float*)&triangle[0].v0,(float*)&triangle[1].v0,(float*)&triangle[2].v0,(float*)&triangle[3].v0);
      const mic_f e1 = gather_4f_zlc(and_mask,(float*)&triangle[0].e1,(float*)&triangle[1].e1,(float*)&triangle[2].e1,(float*)&triangle[3].e1);
      const mic_f e2 = gather_4f_zlc(and_mask,(float*)&triangle[0].e2,(float*)&triangle[1].e2,(float*)&triangle[2].e2,(float*)&triangle[3].e2);
      const mic_f Ng = gather_4f_zlc(and_mask,(float*)&triangle[0].Ng,(float*)&triangle[1].Ng,(float*)&triangle[2].Ng,(float*)&triangle[3].Ng);
      //const mic_f Ng = lcross_xyz(e1,e2);
      
      /* calculate determinant */
      const mic_f C = p0 - org_xyz;
      const mic_f R = lcross_xyz(dir_xyz,C);
      const mic_f det = ldot3_xyz(dir_xyz,Ng);
      const mic_f rcpDet = rcp(det);
      mic_m valid = nz(0x1111,det);
            
      /* perform edge tests */
      const mic_f u = ldot3_xyz(e2,R)*rcpDet;
      const mic_f v = ldot3_xyz(e1,R)*rcpDet;
      valid = gez(valid,u);
      valid = gez(valid,v);
      valid = le (valid,u+v,mic_f(one));
      if (likely(none(valid))) continue;

      /* perform depth test */
      const mic_f t = ldot3_xyz(C,Ng)*rcpDet;
      valid = lt(valid,min_dist_xyz,t);
      valid = gt(valid,max_dist_xyz,t);
      if (likely(none(valid))) continue;

      prefetch<PFHINT_L1EX>(&ray.tfar);
      prefetch<PFHINT_L1EX>(&ray.id0);
      prefetch<PFHINT_L1EX>(&ray.id1);
      prefetch<PFHINT_L1EX>(&ray.u);
      prefetch<PFHINT_L1EX>(&ray.v);
      prefetch<PFHINT_L1EX>(&ray.Ng.x);
      prefetch<PFHINT_L1EX>(&ray.Ng.y);
      prefetch<PFHINT_L1EX>(&ray.Ng.z);

      /* update hit information */
      max_dist_xyz = sel(valid,t,max_dist_xyz);
      const mic_f min_dist = set_min16(max_dist_xyz);
      valid = eq(valid,min_dist,max_dist_xyz);
      max_dist_xyz = min_dist;

      const unsigned long i = bsf64(valid);
      valid = andn(valid,valid-1);
      assert(countbits(valid) == 1);

      const mic_f gnormalx = mic_f(Ng[i+0]);
      const mic_f gnormaly = mic_f(Ng[i+1]);
      const mic_f gnormalz = mic_f(Ng[i+2]);
      
      compactustore16f_low(valid,&ray.u[idx],u);
      compactustore16f_low(valid,&ray.v[idx],v);
      compactustore16f_low(valid,&ray.tfar[idx],t);

      compactustore16f_low(valid,&ray.Ng.x[idx],gnormalx);
      compactustore16f_low(valid,&ray.Ng.y[idx],gnormaly);
      compactustore16f_low(valid,&ray.Ng.z[idx],gnormalz);
      
      const Triangle1* tri_ptr = (Triangle1*)((char*)triangle+16*i);
      ray.id0[idx] = tri_ptr->e1.a;
      ray.id1[idx] = tri_ptr->e2.a;
      
      /* compact the stack */
      if (likely(sindex >= 2))
      {
        if (likely(sindex < 16))
        {
          const unsigned m_num_stack = mic_i::shift1(sindex) - 1;
          const mic_m m_num_stack_low  = toMask(m_num_stack);
          const mic_f snear_low  = upconv16f(stack_near + 0);
          const mic_i snode_low  = upconv16i((int*)stack_node + 0);
          const mic_m m_stack_compact_low  = le(m_num_stack_low,snear_low,max_dist_xyz) | (mic_m)1;
          compactustore16f_low(m_stack_compact_low,stack_near + 0,snear_low);
          compactustore16i_low(m_stack_compact_low,(int*)stack_node + 0,snode_low);
          sindex = countbits(m_stack_compact_low);
          assert(sindex < 16);
        }
        else if (likely(sindex < 32))
        {
          const mic_m m_num_stack_high = toMask(mic_i::shift1(sindex-16) - 1); 
          const mic_f snear_low  = upconv16f(stack_near + 0);
          const mic_f snear_high = upconv16f(stack_near + 16);
          const mic_i snode_low  = upconv16i((int*)stack_node + 0);
          const mic_i snode_high = upconv16i((int*)stack_node + 16);
          const mic_m m_stack_compact_low  = le(snear_low,max_dist_xyz) | (mic_m)1;
          const mic_m m_stack_compact_high = le(m_num_stack_high,snear_high,max_dist_xyz);
          compactustore16f(m_stack_compact_low,      stack_near + 0,snear_low);
          compactustore16i(m_stack_compact_low,(int*)stack_node + 0,snode_low);
          compactustore16f(m_stack_compact_high,      stack_near + countbits(m_stack_compact_low),snear_high);
          compactustore16i(m_stack_compact_high,(int*)stack_node + countbits(m_stack_compact_low),snode_high);
          sindex = countbits(m_stack_compact_low) + countbits(m_stack_compact_high);
          assert ((unsigned)m_num_stack_high == (((unsigned)mic_i::shift1(sindex) - 1) >> 16));
          assert(sindex < 32);
        }
        else if (likely(sindex < 48))
        {
          const mic_m m_num_stack_32 = toMask(mic_i::shift1(sindex-32) - 1); 
          
          const mic_f snear_0  = upconv16f(stack_near + 0);
          const mic_f snear_16 = upconv16f(stack_near + 16);
          const mic_f snear_32 = upconv16f(stack_near + 32);
          const mic_i snode_0  = upconv16i((int*)stack_node + 0);
          const mic_i snode_16 = upconv16i((int*)stack_node + 16);
          const mic_i snode_32 = upconv16i((int*)stack_node + 32);
          const mic_m m_stack_compact_0  = le(               snear_0 ,max_dist_xyz) | (mic_m)1;
          const mic_m m_stack_compact_16 = le(               snear_16,max_dist_xyz);
          const mic_m m_stack_compact_32 = le(m_num_stack_32,snear_32,max_dist_xyz);
          
          sindex = 0;
          compactustore16f(m_stack_compact_0,      stack_near + sindex,snear_0);
          compactustore16i(m_stack_compact_0,(int*)stack_node + sindex,snode_0);
          sindex += countbits(m_stack_compact_0);
          compactustore16f(m_stack_compact_16,      stack_near + sindex,snear_16);
          compactustore16i(m_stack_compact_16,(int*)stack_node + sindex,snode_16);
          sindex += countbits(m_stack_compact_16);
          compactustore16f(m_stack_compact_32,      stack_near + sindex,snear_32);
          compactustore16i(m_stack_compact_32,(int*)stack_node + sindex,snode_32);
          sindex += countbits(m_stack_compact_32);
        }
      }
    }
  }

    static __forceinline bool occluded1(const BVH4AOS* bvh, 
                                 const Node*     const nodes,
                                 const Triangle1 *const triangles,
                                 unsigned rootNode, long idx, 
                                 const mic_f org_xyz, const mic_f dir_xyz, 
                                 const mic_f rdir_xyz, const mic_f org_rdir_xyz, 
                                 const mic_f min_dist_xyz, const mic_f max_dist_xyz)
  {
    /* allocate stack and push root node */
    MIC_ALIGN unsigned stack_node[3*BVH4AOS::maxDepth+1];
    prefetch<PFHINT_L1EX>(&stack_node[0]);      
    prefetch<PFHINT_L1EX>(&stack_node[16]);      
    stack_node[0] = BVH4AOS::invalid_node;
    stack_node[1] = rootNode;
    unsigned sindex = 2;
    const unsigned leaf_mask = BVH4AOS_LEAF_MASK;

    while (1)
    {
      unsigned cur = stack_node[sindex-1];
      sindex--;
      
      while (1) 
      {
        if (unlikely(cur & leaf_mask)) break;

        const Node* __restrict__ const node = (Node*)((char*)nodes + (unsigned long)cur);
        const float* __restrict const plower = (float*)node->lower;
        const float* __restrict const pupper = (float*)node->upper;
        prefetch<PFHINT_L1>((char*)nodes + (unsigned long)cur);
        prefetch<PFHINT_L1>((char*)nodes + (unsigned long)cur + 64);
        
        /* intersect single ray with 4 bounding boxes */
        const mic_f tLowerXYZ = upconv16f(plower) * rdir_xyz - org_rdir_xyz;
        const mic_f tUpperXYZ = upconv16f(pupper) * rdir_xyz - org_rdir_xyz;
        cur = stack_node[sindex-1]; // early pop of next node
        const mic_f tLower = mask_min(0x7777,min_dist_xyz,tLowerXYZ,tUpperXYZ);
        const mic_f tUpper = mask_max(0x7777,max_dist_xyz,tLowerXYZ,tUpperXYZ);
        prefetch<PFHINT_L2>((char*)nodes + (unsigned long)cur);
        prefetch<PFHINT_L2>((char*)nodes + (unsigned long)cur + 64);
        sindex--;
        const mic_f tNear = set_max4(tLower);
        const mic_f tFar  = set_min4(tUpper);  
        const mic_m hitm = le(0x8888,tNear,tFar);
        
        /* if no child is hit, continue with early popped child */
        if (unlikely(ez_mask(hitm) != 0)) continue;

        const unsigned long hiti = toInt(hitm);
        const unsigned long pos_first = bsf64(hiti);
        const unsigned long num_hits = countbits(hiti); 
        sindex++; // undo early pop
        assert(sindex < 3*BVH4AOS::maxDepth+1);

        /* if a single child is hit, continue with that child */
        cur = ((unsigned*)plower)[pos_first];
        if (likely(num_hits == 1)) continue;
        
         /* if two children are hit, push in correct order */
        const unsigned long pos_second = bsf64(pos_first,hiti); 
        if (likely(num_hits == 2))
        {
          const unsigned dist_first = ((unsigned*)&tNear)[pos_first];
          const unsigned dist_second = ((unsigned*)&tNear)[pos_second];
          const unsigned node_second = ((unsigned*)plower)[pos_second];
          
          if (dist_first <= dist_second)
          {
            stack_node[sindex++] = node_second;
            assert(sindex < 3*BVH4AOS::maxDepth+1);
            continue;
          }
          else
          {
            stack_node[sindex++] = cur;
            cur = node_second;
            assert(sindex < 3*BVH4AOS::maxDepth+1);
            continue;
          }
        }

        /* continue with closest child and push all others */
        const mic_f min_dist = set_min_lanes(sel(hitm,tNear,mic_f::inf()));
        const unsigned old_sindex = sindex;
        sindex += countbits(hiti) - 1;
        assert(sindex < 3*BVH4AOS::maxDepth+1);

        const mic_m closest_child = eq(hitm,min_dist,tNear);
        const unsigned long closest_child_pos = bsf64(closest_child);
        const mic_m m_pos = andn(hitm,andn(closest_child,closest_child-1));
        const mic_i plower_node = load16i((int*)plower);
        cur = ((unsigned*)plower)[closest_child_pos];
        compactustore16i(m_pos,&stack_node[old_sindex],plower_node);
      }
      
      /* return if stack is empty */
      if (unlikely(cur == BVH4AOS_INVALID_NODE)) 
        return false;
      
      /* decode leaf node */
      const unsigned itemOffset = BVH4AOS::itemOffset(cur);
      const unsigned items      = BVH4AOS::itemCount (cur);
      const Triangle1* triangle = (Triangle1*)((char*)triangles + itemOffset);

      /* intersect triangles */
      prefetch<PFHINT_L1>(triangle + 3);
      prefetch<PFHINT_L1>(triangle + 2);
      prefetch<PFHINT_L1>(triangle + 1);
      prefetch<PFHINT_L1>(triangle + 0);
      assert(items > 0);
      assert(items <= 4);
      
      /* calculate edges, geometry normal, and determinant */
      const mic_i and_mask = mic_i::zlc4();
      const mic_f p0 = gather_4f_zlc(and_mask,(float*)&triangle[0].v0,(float*)&triangle[1].v0,(float*)&triangle[2].v0,(float*)&triangle[3].v0);
      const mic_f e1 = gather_4f_zlc(and_mask,(float*)&triangle[0].e1,(float*)&triangle[1].e1,(float*)&triangle[2].e1,(float*)&triangle[3].e1);
      const mic_f e2 = gather_4f_zlc(and_mask,(float*)&triangle[0].e2,(float*)&triangle[1].e2,(float*)&triangle[2].e2,(float*)&triangle[3].e2);
      const mic_f Ng = gather_4f_zlc(and_mask,(float*)&triangle[0].Ng,(float*)&triangle[1].Ng,(float*)&triangle[2].Ng,(float*)&triangle[3].Ng);
      //const mic_f Ng = lcross_xyz(e1,e2);
      
      /* calculate determinant */
      const mic_f C = p0 - org_xyz;
      const mic_f R = lcross_xyz(dir_xyz,C);
      const mic_f det = ldot3_xyz(dir_xyz,Ng);
      const mic_f rcpDet = rcp(det);
      mic_m valid = nz(0x1111,det);
            
      /* perform edge tests */
      const mic_f u = ldot3_xyz(e2,R)*rcpDet;
      const mic_f v = ldot3_xyz(e1,R)*rcpDet;
      valid = gez(valid,u);
      valid = gez(valid,v);
      valid = le (valid,u+v,mic_f(one));
      if (likely(none(valid))) continue;

      /* perform depth test */
      const mic_f t = ldot3_xyz(C,Ng)*rcpDet;
      valid = lt(valid,min_dist_xyz,t);
      valid = gt(valid,max_dist_xyz,t);
      if (any(valid)) return true;
    }
    return false;
  }

  private:
    const BVH4AOS* bvh;
  };
}

#endif
