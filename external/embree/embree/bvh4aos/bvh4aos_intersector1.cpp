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

#include "bvh4aos_intersector1.h"
#include "../geometry/triangles.h"

namespace embree
{
  template<typename TriangleIntersector>
  void BVH4AOSIntersector1<TriangleIntersector>::intersect(const BVH4AOSIntersector1* This, Ray& ray)
  {
    /* pointers to node array and triangle array */
    const BVH4AOS* bvh = This->bvh;
    const Node*     const __restrict__ nodes     = (const Node*    ) bvh->nodePtr();
    const Triangle *const __restrict__ triangles = (const Triangle*) bvh->triPtr();
   
    /* allocate stack and push root node */
    MIC_ALIGN unsigned stack_node[3*BVH4AOS::maxDepth+1];
    MIC_ALIGN float    stack_near[3*BVH4AOS::maxDepth+1];
    prefetch<PFHINT_L1EX>(&stack_node[0]);      
    prefetch<PFHINT_L1EX>(&stack_node[16]);      
    prefetch<PFHINT_L1EX>(&stack_near[0]);      
    prefetch<PFHINT_L1EX>(&stack_near[16]);  
    stack_node[0] = BVH4AOS::invalid_node;
    stack_node[1] = bvh->root;
    *(mic_f*)stack_near = mic_f::inf();
    unsigned sindex = 2;
      
    /* load ray into registers */
    const Vector3f ray_rdir  = rcp_safe(ray.dir);
    const mic_f org_xyz      = loadAOS(ray.org.x,ray.org.y,ray.org.z);
    const mic_f dir_xyz      = loadAOS(ray.dir.x,ray.dir.y,ray.dir.z);
    const mic_f rdir_xyz     = loadAOS(ray_rdir.x,ray_rdir.y,ray_rdir.z);
    const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
    mic_f       min_dist_xyz = upconv1f(&ray.tnear);
    mic_f       max_dist_xyz = upconv1f(&ray.tfar);
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
      const unsigned items = BVH4AOS::itemCount(cur);
      const Triangle* triangle = (Triangle*)((char*)triangles + itemOffset);
      
      /* intersect triangles */
      TriangleIntersector::intersect(ray, triangle, items, bvh->vertices);
      max_dist_xyz = upconv1f(&ray.tfar);

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
      }
    }	      
  }

  template<typename TriangleIntersector>
  bool BVH4AOSIntersector1<TriangleIntersector>::occluded(const BVH4AOSIntersector1* This, Ray& ray)
  {
    /* pointers to node array and triangle array */
    const BVH4AOS* bvh = This->bvh;
    const Node*     const __restrict__ nodes     = (const Node*    ) bvh->nodePtr();
    const Triangle *const __restrict__ triangles = (const Triangle*) bvh->triPtr();
   
    /* allocate stack and push root node */
    MIC_ALIGN unsigned stack_node[3*BVH4AOS::maxDepth+1];
    prefetch<PFHINT_L1EX>(&stack_node[0]);      
    prefetch<PFHINT_L1EX>(&stack_node[16]);      
    stack_node[0] = BVH4AOS::invalid_node;
    stack_node[1] = bvh->root;
    unsigned sindex = 2;
      
    /* load ray into registers */
    const Vector3f ray_rdir     = rcp_safe(ray.dir);
    const mic_f org_xyz      = loadAOS(ray.org.x,ray.org.y,ray.org.z);
    const mic_f dir_xyz      = loadAOS(ray.dir.x,ray.dir.y,ray.dir.z);
    const mic_f rdir_xyz     = loadAOS(ray_rdir.x,ray_rdir.y,ray_rdir.z);
    const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
    const mic_f min_dist_xyz = upconv1f(&ray.tnear);
    const mic_f max_dist_xyz = upconv1f(&ray.tfar);
    const unsigned leaf_mask = BVH4AOS_LEAF_MASK;

    while (1)
    {
      NodeRef cur = stack_node[--sindex];
      
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
      const unsigned items = BVH4AOS::itemCount(cur);
      const Triangle* triangle = (Triangle*)((char*)triangles + itemOffset);
      
      /* intersect triangles */
      if (TriangleIntersector::occluded(ray, triangle, items, bvh->vertices))
        return true;
    }
    return false;
  }

  void BVH4AOSIntersector1Register () 
  {
    //TriangleMesh::intersectors1.add("bvh4aos","triangle1i","fast","moeller" ,false,BVH4AOSIntersector1<Triangle1iIntersector1<Intersector1MoellerTrumbore> >::create);
    //TriangleMesh::intersectors1.add("bvh4aos","triangle1i","fast","pluecker",true ,BVH4AOSIntersector1<Triangle1iIntersector1<Intersector1Pluecker> >::create);
    //TriangleMesh::intersectors1.add("bvh4aos","triangle1v","fast","moeller" ,false,BVH4AOSIntersector1<Triangle1vIntersector1<Intersector1MoellerTrumbore> >::create);
    //TriangleMesh::intersectors1.add("bvh4aos","triangle1v","fast","pluecker",true ,BVH4AOSIntersector1<Triangle1vIntersector1<Intersector1Pluecker> >::create);
    //TriangleMesh::intersectors1.add("bvh4aos","triangle4i","fast","moeller" ,false,BVH4AOSIntersector1<Triangle4iIntersector1<Intersector1MoellerTrumbore> >::create);
    //TriangleMesh::intersectors1.add("bvh4aos","triangle4i","fast","pluecker",true ,BVH4AOSIntersector1<Triangle4iIntersector1<Intersector1Pluecker> >::create);
    //TriangleMesh::intersectors1.add("bvh4aos","triangle4v","fast","moeller" ,false,BVH4AOSIntersector1<Triangle4vIntersector1<Intersector1MoellerTrumbore> >::create);
    //TriangleMesh::intersectors1.add("bvh4aos","triangle4v","fast","pluecker",true ,BVH4AOSIntersector1<Triangle4vIntersector1<Intersector1Pluecker> >::create);
    //TriangleMesh::intersectors1.add("bvh4aos","triangle4aos","fast","moeller",true ,BVH4AOSIntersector1<Triangle4AOSIntersector1MoellerTrumboreMIC>::create);
    TriangleMesh::intersectors1.add("bvh4aos","triangle1","fast","moeller",true ,BVH4AOSIntersector1<Triangle1Intersector1MoellerTrumboreMIC>::create);
    //TriangleMesh::intersectors1.add("bvh4aos","triangle4" ,"fast","moeller" ,true ,BVH4AOSIntersector1<Triangle4Intersector1MoellerTrumbore>::create);
    //TriangleMesh::intersectors1.add("bvh4aos","triangle8" ,"fast","moeller" ,true ,BVH4AOSIntersector1<Triangle8Intersector1MoellerTrumbore>::create);
    TriangleMesh::intersectors1.setAccelDefaultTraverser("bvh4aos","fast");

    VirtualScene::intersectors1.add("bvh4aos","virtual","fast","virtual" ,true ,BVH4AOSIntersector1<VirtualObjectIntersector1>::create);
    VirtualScene::intersectors1.setAccelDefaultTraverser("bvh4aos","fast");
  }
}

