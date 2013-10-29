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

#include "bvh4aos_triangle1_intersector16_hybrid.h"
#include "bvh4aos_triangle1_intersector1.h"
#include "geometry/triangles.h"

#define SIMD_UTIL_SWITCH_THRESHOLD 8

namespace embree
{
  void BVH4AOSTriangle1Intersector16Hybrid::intersect(const BVH4AOSTriangle1Intersector16Hybrid* This, Ray16& ray, const __mmask valid_i)
  {
    /* pointers to node array and triangle array */
    const mic_m valid = valid_i;
    const BVH4AOS* bvh = This->bvh;
    const Node*     const __restrict__ nodes     = (const Node*    ) bvh->nodePtr();
    const Triangle* const __restrict__ triangles = (const Triangle*) bvh->triPtr();

    /* load ray into registers */
    mic_f min_distance = sel(valid,ray.tnear,pos_inf);
    mic_f max_distance = sel(valid,ray.tfar ,neg_inf);
    
    /* allocate stack and push root node */
    MIC_ALIGN unsigned stack_node[3*BVH4AOS::maxDepth+1];
    MIC_ALIGN mic_f    stack_near[3*BVH4AOS::maxDepth+1];
    prefetch<PFHINT_L1EX>(stack_node);      
    prefetch<PFHINT_L1EX>(stack_near);      
    stack_node[0] = BVH4AOS::invalid_node;
    stack_near[0] = mic_f::inf();
    stack_node[1] = bvh->root;
    stack_near[1] = min_distance;
    unsigned* sptr_node = stack_node+2;
    mic_f   * sptr_near = stack_near+2;

    /* load ray into registers */
    const mic_f orgx = ray.org.x;
    const mic_f orgy = ray.org.y;
    const mic_f orgz = ray.org.z; 
    const mic_f dirx = ray.dir.x;
    const mic_f diry = ray.dir.y;
    const mic_f dirz = ray.dir.z; 
    const mic_f rdirx = rcp_safe(ray.dir.x);
    const mic_f rdiry = rcp_safe(ray.dir.y);
    const mic_f rdirz = rcp_safe(ray.dir.z);
    const mic_f org_rdirx = orgx * rdirx;
    const mic_f org_rdiry = orgy * rdiry;
    const mic_f org_rdirz = orgz * rdirz;
    const unsigned leaf_mask = BVH4AOS_LEAF_MASK;

    while (1)
    {
      /* pop next node from stack */
      mic_f minDist = upconv16f((const float*)(sptr_near-1));
      const mic_m m_dist = gt(max_distance,upconv16f((const float*)(sptr_near-1)));
      unsigned cur = *(sptr_node-1);
      sptr_node--;
      sptr_near--;
      if (unlikely(cur == BVH4AOS_INVALID_NODE)) 
        break;

      /* cull node if behind closest hit point */
      if (unlikely(none(m_dist))) continue;

      /* trace single rays */
      if (unlikely(countbits(m_dist) <= SIMD_UTIL_SWITCH_THRESHOLD)) 
      { 
        store16f(MIC_M_ALL,&ray.tfar,max_distance);
        long rayIndex = -1;
        const unsigned int m_index = toInt(m_dist);
        while((rayIndex = bsf64(rayIndex,m_index)) != MIC_NO_BIT_SET_64)	    
        {	    
          // === TODO: precompute SOAtoAOS transformation, load with 4x broadcast
          const mic_f org_xyz      = SOAtoAOS_4f(rayIndex,ray.org.x,ray.org.y,ray.org.z);
          const mic_f dir_xyz      = SOAtoAOS_4f(rayIndex,dirx,diry,dirz);
          const mic_f rdir_xyz     = SOAtoAOS_4f(rayIndex,rdirx,rdiry,rdirz);
          const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
          const mic_f min_dist_xyz = upconv1f(&min_distance[rayIndex]); 
          const mic_f max_dist_xyz = upconv1f(&max_distance[rayIndex]); 
          BVH4AOSTriangle1Intersector1::intersect1(bvh,nodes,triangles,cur,rayIndex,org_xyz,dir_xyz,rdir_xyz,org_rdir_xyz,min_dist_xyz,max_dist_xyz,ray);
        }
        max_distance = ray.tfar;
        continue; 
      }

      while (1)
      {
        /* test if this is a leaf node */
        if (unlikely(cur & leaf_mask)) 
          break;
        
        /* pop of next node */
        sptr_node--;
        sptr_near--;
        const Node* const __restrict__ node = (Node*)((char*)nodes + (unsigned long)cur);
        prefetch<PFHINT_L1>((mic_f*)node + 1); // depth first order prefetch		       
        cur = *sptr_node;
        minDist = *sptr_near;
                
#pragma unroll(4)
        for (unsigned i=0;i<4;i++)
	{
          const float* const __restrict__ plower = (float*) &node->lower[i];
          const float* const __restrict__ pupper = (float*) &node->upper[i];
          
          const mic_f tLowerX = plower[0] * rdirx - org_rdirx;
          const mic_f tLowerY = plower[1] * rdiry - org_rdiry;
          const mic_f tLowerZ = plower[2] * rdirz - org_rdirz;
          const mic_f tUpperX = pupper[0] * rdirx - org_rdirx;
          const mic_f tUpperY = pupper[1] * rdiry - org_rdiry;
          const mic_f tUpperZ = pupper[2] * rdirz - org_rdirz;
          
          const mic_f tNear = _max(_max(_min(tLowerX, tUpperX), _min(tLowerY, tUpperY)), _min(tLowerZ, tUpperZ));
          const mic_f tFar  = _min(_min(_max(tLowerX, tUpperX), _max(tLowerY, tUpperY)), _max(tLowerZ, tUpperZ));      
          const mic_m mhit   = le(_max(min_distance,tNear), _min(tFar,max_distance));
          const mic_f childDist = sel(mhit,tNear,inf);
          const mic_m closer = lt(childDist,minDist);
          
          /* if we hit the child we choose to continue with that child if it 
             is closer than the current next child, or we push it onto the stack */
          if (likely(any(mhit)))
          {
            const unsigned child = node->lower[i].child;
            sptr_node++;
            sptr_near++;
            
            /* push cur node onto stack and continue with hit child */
            if (any(closer)) {
              *(sptr_node-1) = cur;
              *(sptr_near-1) = minDist; 
              minDist = childDist;
              cur = child;
            } 

            /* push hit child onto stack*/
            else {
              *(sptr_node-1) = child;
              *(sptr_near-1) = childDist; 
            }          
            assert(sptr_node - stack_node < 3*BVH4AOS::maxDepth+1);
          }	      
        }
      }
      
      /* return if stack is empty */
      if (unlikely(cur == BVH4AOS_INVALID_NODE)) 
        break;
      
      /* decode leaf node */
      const unsigned itemOffset = BVH4AOS::itemOffset(cur);
      const unsigned items      = BVH4AOS::itemCount (cur);
      const Triangle1* triangle = (Triangle1*)((char*)triangles + itemOffset);
      assert(items > 0);
      assert(items <= 4);
      
      /* intersect triangles */
      mic3f org(orgx,orgy,orgz);
      mic3f dir(dirx,diry,dirz);
      prefetch<PFHINT_L1>((mic_f*)triangle +  0); 
      prefetch<PFHINT_L2>((mic_f*)triangle +  1); 
      prefetch<PFHINT_L2>((mic_f*)triangle +  2); 
      prefetch<PFHINT_L2>((mic_f*)triangle +  3); 

      for (unsigned int i=0; i<items; i++, triangle++) 
      {
        prefetch<PFHINT_L1>((mic_f*)triangle+1); 
        
        /* calculate edges and geometry normal */
        const mic_f _p0 = upconv4f((float*)&triangle->v0);
        const mic_f _e1 = upconv4f((float*)&triangle->e1);
        const mic_f _e2 = upconv4f((float*)&triangle->e2);
        const mic_f _Ng = upconv4f((float*)&triangle->Ng);
        //const mic_f _Ng = lcross_xyz(_e1,_e2);
        
        /* change data layout */
        const mic3f p0(swAAAA(_p0),swBBBB(_p0),swCCCC(_p0));
        const mic3f e1(swAAAA(_e1),swBBBB(_e1),swCCCC(_e1));
        const mic3f e2(swAAAA(_e2),swBBBB(_e2),swCCCC(_e2));
        const mic3f Ng(swAAAA(_Ng),swBBBB(_Ng),swCCCC(_Ng));
        
        /* calculate determinant */
        const mic3f C = p0 - org;
        const mic3f R = cross(dir,C);
        const mic_f det = dot(dir,Ng);
        const mic_f rcpDet = rcp(det);
        mic_m valid1 = nz(valid,det);
        
        /* perform edge tests */
        const mic_f u = dot(e2,R)*rcpDet;
        const mic_f v = dot(e1,R)*rcpDet;
        valid1 = gez(valid1,u);
        valid1 = gez(valid1,v);
        valid1 = le (valid1,u+v,mic_f(one));
        if (likely(none(valid1))) continue;
        
        /* perform depth test */
        const mic_f t = dot(C,Ng)*rcpDet;
        valid1 = lt(valid1,min_distance,t);
        valid1 = gt(valid1,max_distance,t);
        const mic_i id0 = upconv1i((const int *)&triangle->e1.a);
        const mic_i id1 = upconv1i((const int *)&triangle->e2.a);
        if (likely(none(valid1))) continue;
        
        /* update hit information */
        prefetch<PFHINT_L1EX>(&ray.u);
        prefetch<PFHINT_L1EX>(&ray.v);
        prefetch<PFHINT_L1EX>(&ray.tfar);
        prefetch<PFHINT_L1EX>(&ray.id0);
        prefetch<PFHINT_L1EX>(&ray.id1);
        max_distance = sel(valid1,t,max_distance);
        store16f(valid1,&ray.Ng.x,Ng.x);
        store16f(valid1,&ray.Ng.y,Ng.y);
        store16f(valid1,&ray.Ng.z,Ng.z);
        store16f(valid1,&ray.u,u);
        store16f(valid1,&ray.v,v);
        store16f(valid1,&ray.tfar,t);
        store16i(valid1,&ray.id0,id0);
        store16i(valid1,&ray.id1,id1);
      }
    }
  }

  __mmask BVH4AOSTriangle1Intersector16Hybrid::occluded(const BVH4AOSTriangle1Intersector16Hybrid* This, Ray16& ray, const __mmask valid_i)
  {
    /* pointers to node array and triangle array */
    const mic_m valid = valid_i;
    const BVH4AOS* bvh = This->bvh;
    const Node*     const __restrict__ nodes     = (const Node*    ) bvh->nodePtr();
    const Triangle* const __restrict__ triangles = (const Triangle*) bvh->triPtr();

    /* load ray into registers */
    mic_f min_distance = sel(valid,ray.tnear,pos_inf);
    mic_f max_distance = sel(valid,ray.tfar ,neg_inf);
    
    /* allocate stack and push root node */
    MIC_ALIGN unsigned stack_node[3*BVH4AOS::maxDepth+1];
    MIC_ALIGN mic_f    stack_near[3*BVH4AOS::maxDepth+1];
    prefetch<PFHINT_L1EX>(stack_node);      
    prefetch<PFHINT_L1EX>(stack_near);      
    stack_node[0] = BVH4AOS::invalid_node;
    stack_near[0] = mic_f::inf();
    stack_node[1] = bvh->root;
    stack_near[1] = min_distance;
    unsigned* sptr_node = stack_node+2;
    mic_f   * sptr_near = stack_near+2;

    /* load ray into registers */
    const mic_f orgx = ray.org.x;
    const mic_f orgy = ray.org.y;
    const mic_f orgz = ray.org.z; 
    const mic_f dirx = ray.dir.x;
    const mic_f diry = ray.dir.y;
    const mic_f dirz = ray.dir.z; 
    const mic_f rdirx = rcp_safe(ray.dir.x);
    const mic_f rdiry = rcp_safe(ray.dir.y);
    const mic_f rdirz = rcp_safe(ray.dir.z);
    const mic_f org_rdirx = orgx * rdirx;
    const mic_f org_rdiry = orgy * rdiry;
    const mic_f org_rdirz = orgz * rdirz;
    const unsigned leaf_mask = BVH4AOS_LEAF_MASK;
    mic_m m_not_occluded = valid;

    while (1)
    {
      /* pop next node from stack */
      mic_f minDist = upconv16f((const float*)(sptr_near-1));
      unsigned cur = *(sptr_node-1);
      sptr_node--;
      sptr_near--;
      if (unlikely(cur == BVH4AOS_INVALID_NODE)) 
        break;

      /* cull node if behind closest hit point */
      const mic_m m_dist = gt(m_not_occluded,max_distance,minDist);
      if (unlikely(none(m_dist))) continue;

      /* trace single rays */
      if (unlikely(countbits(m_dist) <= SIMD_UTIL_SWITCH_THRESHOLD)) 
      { 
        mic_i not_occluded = sel(m_not_occluded,mic_i::minus_one(),mic_i::zero());
        long rayIndex = -1;
        const unsigned int m_index = toInt(m_dist);
        while((rayIndex = bsf64(rayIndex,m_index)) != MIC_NO_BIT_SET_64)	    
        {	    
          // === TODO: precompute SOAtoAOS transformation, load with 4x broadcast
          const mic_f org_xyz      = SOAtoAOS_4f(rayIndex,ray.org.x,ray.org.y,ray.org.z);
          const mic_f dir_xyz      = SOAtoAOS_4f(rayIndex,dirx,diry,dirz);
          const mic_f rdir_xyz     = SOAtoAOS_4f(rayIndex,rdirx,rdiry,rdirz);
          const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
          const mic_f min_dist_xyz = upconv1f(&min_distance[rayIndex]); 
          const mic_f max_dist_xyz = upconv1f(&max_distance[rayIndex]); 
          if (BVH4AOSTriangle1Intersector1::occluded1(bvh,nodes,triangles,cur,rayIndex,org_xyz,dir_xyz,rdir_xyz,org_rdir_xyz,min_dist_xyz,max_dist_xyz))
            not_occluded[rayIndex] = 0;
        }
        m_not_occluded = eq(not_occluded,mic_i::minus_one());
        if (unlikely(m_not_occluded == 0)) {
          return valid;
        }
        continue; 
      }

      while (1)
      {
        /* test if this is a leaf node */
        if (unlikely(cur & leaf_mask)) 
          break;
        
        /* pop of next node */
        sptr_node--;
        sptr_near--;
        const Node* const node = (Node*)((char*)nodes + (unsigned long)cur);
        prefetch<PFHINT_L1>((mic_f*)node + 1); // depth first order prefetch		
        cur = *sptr_node;
        minDist = *sptr_near;
                
#pragma unroll(4)
        for (unsigned i=0;i<4;i++)
	{
          const float* const plower = (float*) &node->lower[i];
          const float* const pupper = (float*) &node->upper[i];
          
          const mic_f tLowerX = plower[0] * rdirx - org_rdirx;
          const mic_f tLowerY = plower[1] * rdiry - org_rdiry;
          const mic_f tLowerZ = plower[2] * rdirz - org_rdirz;
          const mic_f tUpperX = pupper[0] * rdirx - org_rdirx;
          const mic_f tUpperY = pupper[1] * rdiry - org_rdiry;
          const mic_f tUpperZ = pupper[2] * rdirz - org_rdirz;
          
          const mic_f tNear = _max(_max(_min(tLowerX, tUpperX), _min(tLowerY, tUpperY)), _min(tLowerZ, tUpperZ));
          const mic_f tFar  = _min(_min(_max(tLowerX, tUpperX), _max(tLowerY, tUpperY)), _max(tLowerZ, tUpperZ));      
          const mic_m mhit   = le(m_not_occluded,_max(min_distance,tNear), _min(tFar,max_distance));
          const mic_f childDist = sel(mhit,tNear,inf);
          const mic_m closer = lt(m_not_occluded,childDist,minDist);
          
          /* if we hit the child we choose to continue with that child if it 
             is closer than the current next child, or we push it onto the stack */
          if (likely(any(mhit)))
          {
            const unsigned child = node->lower[i].child;
            sptr_node++;
            sptr_near++;
            
            /* push cur node onto stack and continue with hit child */
            if (any(closer)) {
              *(sptr_node-1) = cur;
              *(sptr_near-1) = minDist; 
              minDist = childDist;
              cur = child;
            } 

            /* push hit child onto stack*/
            else {
              *(sptr_node-1) = child;
              *(sptr_near-1) = childDist; 
            }          
            assert(sptr_node - stack_node < 3*BVH4AOS::maxDepth+1);
          }	      
        }
      }
      
      /* return if stack is empty */
      if (unlikely(cur == BVH4AOS_INVALID_NODE))
        break;
      
      /* decode leaf node */
      const unsigned itemOffset = BVH4AOS::itemOffset(cur);
      const unsigned items      = BVH4AOS::itemCount (cur);
      const Triangle1* triangle = (Triangle1*)((char*)triangles + itemOffset);

      /* intersect triangles */
      mic_m valid_leaf = lt(m_not_occluded,minDist,max_distance);

      assert(items > 0);
      assert(items <= 4);
      prefetch<PFHINT_L1>((mic_f*)triangle +  0); 
      prefetch<PFHINT_L2>((mic_f*)triangle +  1); 
      prefetch<PFHINT_L2>((mic_f*)triangle +  2); 
      prefetch<PFHINT_L2>((mic_f*)triangle +  3); 
      mic3f org(orgx,orgy,orgz);
      mic3f dir(dirx,diry,dirz);

      for (unsigned int i=0; i<items;i++,triangle++) 
      {
        prefetch<PFHINT_L1>((mic_f*)triangle +  1); 
        
        /* calculate edges and geometry normal */
        const mic_f _p0 = upconv4f((float*)&triangle->v0);
        const mic_f _e1 = upconv4f((float*)&triangle->e1);
        const mic_f _e2 = upconv4f((float*)&triangle->e2);
        const mic_f _Ng = upconv4f((float*)&triangle->Ng);
        //const mic_f _Ng = lcross_xyz(_e1,_e2);
        
        /* change data layout */
        const mic3f p0(swAAAA(_p0),swBBBB(_p0),swCCCC(_p0));
        const mic3f e1(swAAAA(_e1),swBBBB(_e1),swCCCC(_e1));
        const mic3f e2(swAAAA(_e2),swBBBB(_e2),swCCCC(_e2));
        const mic3f Ng(swAAAA(_Ng),swBBBB(_Ng),swCCCC(_Ng));
        
        /* calculate determinant */
        const mic3f C = p0 - org;
        const mic3f R = cross(dir,C);
        const mic_f det = dot(dir,Ng);
        const mic_f rcpDet = rcp(det);
        mic_m valid1 = nz(valid_leaf,det);
        
        /* perform edge tests */
        const mic_f u = dot(e2,R)*rcpDet;
        const mic_f v = dot(e1,R)*rcpDet;
        valid1 = gez(valid1,u);
        valid1 = gez(valid1,v);
        valid1 = le (valid1,u+v,mic_f(one));
        if (likely(none(valid1))) continue;
        
        /* perform depth test */
        const mic_f t = dot(C,Ng)*rcpDet;
        valid1 = lt(valid1,min_distance,t);
        valid1 = gt(valid1,max_distance,t);

        m_not_occluded &= ~valid1;
        valid_leaf &= ~valid1;
        if (unlikely(valid_leaf == 0)) break;
      }
      if (unlikely(m_not_occluded == 0)) {
        return valid;
      }
    }
    return valid & ~m_not_occluded;
  }

  void BVH4AOSTriangle1Intersector16HybridRegister () 
  {
    TriangleMesh::intersectors16.add("bvh4aos","triangle1","hybrid","moeller",true ,BVH4AOSTriangle1Intersector16Hybrid::create);
  }
}
