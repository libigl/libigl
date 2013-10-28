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

#include "bvh4aos_triangle1_intersector16_chunk_ref.h"

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

  void BVH4AOSTriangle1Intersector16ChunkRef::intersect(const BVH4AOSTriangle1Intersector16ChunkRef* This, Ray16& ray, const __mmask m_active_i) //(const mic_m m_active, Ray &ray)
  {
    MIC_ALIGN unsigned int stack_node[QUAD_BVH_MAX_STACK_DEPTH*2];
    MIC_ALIGN mic_f stack_near[QUAD_BVH_MAX_STACK_DEPTH*2];
    
    //const mic3f &origin    = ray.org;
    //const mic3f &direction = ray.dir;
    const mic_f min_dist = ray.tnear; //const mic_f &min_dist  = ray.t_min;
    const mic_f max_dist = ray.tfar; //const mic_f max_dist   = ray.t;
    
    const mic_m m_active = m_active_i;
    const BVH4AOS* bvh = This->bvh;
    const Node      *const __restrict__ qbvh  = (const Node*) bvh->nodePtr(); //Renderer::scene->qbvh;
    const Triangle1 *const __restrict__ accel = (const Triangle1*) bvh->triPtr(); //Renderer::scene->accel;
    
    const mic_f orgx = ray.org.x;
    const mic_f orgy = ray.org.y;
    const mic_f orgz = ray.org.z;  
    
    //mic_i triangleID = mic_i::minus_one();
    //mic_i shaderID   = mic_i::minus_one();
    mic_i ID0 = mic_i::minus_one();
    mic_i ID1 = mic_i::minus_one();
    
    mic_f min_distance   = min_dist;
    mic_f max_distance   = sel(m_active,max_dist,mic_f::zero());

    const mic_f direction_x = sel(eqz(ray.dir.x),mic_f::ulp(),ray.dir.x);
    const mic_f direction_y = sel(eqz(ray.dir.y),mic_f::ulp(),ray.dir.y);
    const mic_f direction_z = sel(eqz(ray.dir.z),mic_f::ulp(),ray.dir.z);
    
    const mic_f rdirx = rcp(direction_x);
    const mic_f rdiry = rcp(direction_y);
    const mic_f rdirz = rcp(direction_z);

    prefetch<PFHINT_L1EX>(stack_node);      
    prefetch<PFHINT_L1EX>(stack_near);      
    
    const mic_f org_rdirx = orgx * rdirx;
    const mic_f org_rdiry = orgy * rdiry;
    const mic_f org_rdirz = orgz * rdirz;
    
    const mic_f inf = mic_f::inf();
    const mic_f zero = mic_f::zero();
    
    /* --  dummy node at index 0, removes branch in inner-loop -- */
    stack_node[0] = QBVH_TERMINAL_TOKEN;
    stack_near[0] = inf;

    unsigned int* __restrict__ sptr_node = stack_node + 1;
    mic_f       * __restrict__ sptr_near = stack_near + 1;

    *sptr_node++ = bvh->root; //qbvh[0].lower[0].d;
    *sptr_near++ = min_dist; // mic_f::eps();

    STAT_COLLECTOR(StatCollector::chunk.numChunks++; StatCollector::chunk.numRays += countbits(m_active));

    while (1)
    {
      mic_f minDist = upconv16f((const float*)(sptr_near-1));
      const mic_m m_dist = gt(max_distance,upconv16f((float*)(sptr_near-1)));
      unsigned int curNode = *(sptr_node-1);

      //const mic_m m_dist = lt(minDist,max_distance);

      sptr_node--;
      sptr_near--;

      if (unlikely(curNode == QBVH_TERMINAL_TOKEN)) break;

      STAT_COLLECTOR(StatCollector::chunk.numStackPops++);

      if (unlikely(m_dist == 0)) continue;

      assert(curNode != QBVH_TERMINAL_TOKEN);


      while(1)
      {
        if (unlikely(qbvhLeaf(curNode) != 0)) break;
        
        assert (curNode != QBVH_TERMINAL_TOKEN);
        
        //mic_f minDist = *(sptr_near-1);
        
        sptr_node--;
        sptr_near--;
        
        const Node* __restrict__ const qptr = qbvhChildPtr(qbvh,curNode);
        
        STAT_COLLECTOR(StatCollector::chunk.numTraversalSteps++);
        prefetch<PFHINT_L1>((mic_f*)qptr + 1); // depth first order, prefetch		
        
        STAT_COLLECTOR(unsigned int quadhits=0);
        curNode = *sptr_node;
        minDist = *sptr_near;

        STAT_COLLECTOR(StatCollector::chunk.numActiveRaysTraversal+=countbits(m_dist));
        
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
          
          const mic_m lhit   = le(_max(min_distance,lnearP), _min(lfarP,max_distance));
          
          const mic_f boxDist = sel(lhit,lnearP,inf);
          
          const mic_m m_update_min = lt(boxDist,minDist);
          
          if (likely(nz_mask(lhit) != 0))
          {
            const unsigned int dd = qptr->lower[i].child;
            sptr_node++;
            sptr_near++;
            
            if (nz_mask(m_update_min)  != 0) {
              *(sptr_node-1) = curNode;
              *(sptr_near-1) = minDist; 
              minDist = boxDist;
              curNode = dd;
            } else {
              *(sptr_node-1) = dd;
              *(sptr_near-1) = boxDist; 
            }
            
            assert(sptr_node - stack_node < QUAD_BVH_MAX_STACK_DEPTH);
          }	      
        }
        STAT_COLLECTOR(StatCollector::chunk.histo[quadhits]++);
      }
      
      if (unlikely(curNode == QBVH_TERMINAL_TOKEN)) break;
      assert (curNode != QBVH_TERMINAL_TOKEN);
      
      STAT_COLLECTOR(StatCollector::chunk.numActiveRaysAtLeaf+=countbits(lt(minDist,max_distance)));
      STAT_COLLECTOR(StatCollector::chunk.numLeafIntersections++);
      
      unsigned int itemOffset    = qbvhItemOffset(curNode);
      const Triangle1  *__restrict__ tptr    = (Triangle1*)((char*)accel + itemOffset);
      const unsigned int items = qbvhItems(curNode);
      unsigned int triID = qbvhItemOffsetToID(itemOffset);
      
      assert(items > 0);
      assert(items <= 4);

      prefetch<PFHINT_L1>((mic_f*)tptr +  0); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  1); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  2); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  3); 


      STAT_COLLECTOR(StatCollector::chunk.numActiveRaysIntersection+=items * countbits(lt(minDist,max_distance)));

      for (unsigned int i=0;i<items;i++,tptr++,triID++) 
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
        
        //const mic_f normal = lcross_zxy(e1,e2);
        const mic_f normal = lcross_zxy(e2,e1);
        const mic_f nz = swAAAA(normal);
        const mic_f nx = swBBBB(normal);
        const mic_f ny = swCCCC(normal);
        
        const mic_f nom = madd231(madd231(v0_org_z * nz,v0_org_y, ny),v0_org_x,nx);
        const mic_f den = madd231(madd231(direction_z * nz,direction_y,ny),direction_x,nx);
        
#if defined(BACK_FACE_CULLING)    
        const mic_m m_sign = lz(den);
        if (unlikely(m_sign == 0)) continue;    
#endif
        
        const mic_f odz = msubr231(direction_x*v0_org_y,v0_org_x,direction_y);
        const mic_f odx = msubr231(direction_y*v0_org_z,v0_org_y,direction_z);
        const mic_f ody = msubr231(direction_z*v0_org_x,v0_org_z,direction_x);	  
        
        const mic_f u = madd231(madd231(odz*swCCCC(e2),ody,swBBBB(e2)),odx,swAAAA(e2));
        const mic_f v = madd231(madd231(odz*swCCCC(e1),ody,swBBBB(e1)),odx,swAAAA(e1));
        
#if !defined(BACK_FACE_CULLING)    
        const mic_m den_gt_0 = gz(den); 
        
        const mic_m valid_u = (den_gt_0^ge(u,zero));
        const mic_m valid_v = (den_gt_0^ge(v,zero));
        const mic_m valid_uv = den_gt_0^ge(add_neg(u,v),den); 
        const mic_m m_aperture = valid_u & valid_v & valid_uv;
#else
        const mic_m den_gt_0 = lt(den,zero);
        const mic_m valid_u = ge(den_gt_0,u,zero); 
        const mic_m valid_v = ge(valid_u,v,zero); 
        const mic_m m_aperture = ge(valid_v,add_neg(u,v),den); 
#endif
        
        const mic_f rcp_den = rcp(den);
        if ( unlikely(m_aperture == 0) ) continue;
        
        const mic_i i_triID = mic_i(triID);
        //const mic_i shdID = upconv1i((const int *)&tptr->shaderID);
        const mic_i id0 = upconv1i((const int *)&tptr->e1.a);
        const mic_i id1 = upconv1i((const int *)&tptr->e2.a);
        
        const mic_f t = rcp_den*nom;
        const mic_m m_final  = lt(lt(m_aperture,min_distance,t),t,max_distance);
        
        //prefetch<PFHINT_L1EX>(&ray.gNormal[0]);      
        //prefetch<PFHINT_L1EX>(&ray.gNormal[1]);      
        //prefetch<PFHINT_L1EX>(&ray.gNormal[2]);      
        
        max_distance = sel(m_final,t,max_distance);
        
        //triangleID   = sel(m_final,i_triID,triangleID); 
        //shaderID     = sel(m_final,shdID,shaderID); 
        ID0 = sel(m_final,id0,ID0); 
        ID1 = sel(m_final,id1,ID1); 
        
        prefetch<PFHINT_L1EX>(&ray.u); //prefetch<PFHINT_L1EX>(&ray.u);      
        prefetch<PFHINT_L1EX>(&ray.v); //prefetch<PFHINT_L1EX>(&ray.v);      
        
        if ( unlikely(m_final == 0) ) continue;
        
        prefetch<PFHINT_L1EX>(&ray.id0); //prefetch<PFHINT_L1EX>(&ray.triID);      
        prefetch<PFHINT_L1EX>(&ray.id1); //prefetch<PFHINT_L1EX>(&ray.shaderID);      
        prefetch<PFHINT_L1EX>(&ray.tfar); //prefetch<PFHINT_L1EX>(&ray.t);      
        
        //const mic_f rden = sel(den_gt_0,-rcp_den,rcp_den);
        //const mic_f rden = neg(~den_gt_0,rcp_den,rcp_den);
        //const mic_f rden = rcp_den;
        const mic_f du = abs(u*rcp_den);		      
        const mic_f dv = abs(v*rcp_den);
        
        store16f(m_final,&ray.Ng.x,nx); //store16f(m_final,&ray.gNormal[0],nx);
        store16f(m_final,&ray.Ng.y,ny); //store16f(m_final,&ray.gNormal[1],ny);
        store16f(m_final,&ray.Ng.z,nz); //store16f(m_final,&ray.gNormal[2],nz);
        
        store16f(m_final,&ray.u,du); //store16f(m_final,&ray.u,du);
        store16f(m_final,&ray.v,dv); //store16f(m_final,&ray.v,dv);
      }
    }
    store16f(MIC_M_ALL,&ray.tfar,max_distance); //store16f(MIC_M_ALL,&ray.t,max_distance);
    store16i(MIC_M_ALL,&ray.id0,ID0);        //store16i(MIC_M_ALL,&ray.triID,triangleID); 
    store16i(MIC_M_ALL,&ray.id1,ID1);        //store16i(MIC_M_ALL,&ray.shaderID,shaderID); 
  }

  __mmask BVH4AOSTriangle1Intersector16ChunkRef::occluded(const BVH4AOSTriangle1Intersector16ChunkRef* This, Ray16& ray, const __mmask m_active_i) //(const mic_m m_active, Ray &ray)
  {
    MIC_ALIGN unsigned int stack_node[QUAD_BVH_MAX_STACK_DEPTH*2];
    MIC_ALIGN mic_f stack_near[QUAD_BVH_MAX_STACK_DEPTH*2];
    
    const mic_m m_active = m_active_i;
    const BVH4AOS* bvh = This->bvh;
    const Node      *const __restrict__ qbvh  = (const Node*) bvh->nodePtr(); //Renderer::scene->qbvh;
    const Triangle1 *const __restrict__ accel = (const Triangle1*) bvh->triPtr(); //Renderer::scene->accel;
       
    //const mic3f &origin    = ray.org;
    //const mic3f &direction = ray.dir;
    const mic_f min_dist  = ray.tnear;
    const mic_f max_dist   = ray.tfar;
    
    const mic_f orgx = ray.org.x; //origin.x;
    const mic_f orgy = ray.org.y; //origin.y;
    const mic_f orgz = ray.org.z; //origin.z;  
    
    const mic_f min_distance   = min_dist;
    const mic_f max_distance   = sel(m_active,max_dist,mic_f::zero());
    
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
    const mic_f zero = mic_f::zero();
    
    /* --  dummy node at index 0, removes branch in inner-loop -- */
    stack_node[0] = QBVH_TERMINAL_TOKEN;
    stack_near[0] = inf;
    
    unsigned int* __restrict__ sptr_node = stack_node + 1;
    mic_f       * __restrict__ sptr_near = stack_near + 1;
    
    *sptr_node++ = bvh->root; //qbvh[0].lower[0].d;
    *sptr_near++ = min_dist; // mic_f::eps();
    
    STAT_COLLECTOR(StatCollector::chunk.numChunks++; StatCollector::chunk.numRays += countbits(m_active));
    
    mic_m m_not_occluded = m_active;
    
    while (1)
    {
      mic_f minDist = upconv16f((const float*)(sptr_near-1));
      const mic_m m_dist = gt(m_not_occluded,max_distance,upconv16f((float*)(sptr_near-1)));
      unsigned int curNode = *(sptr_node-1);
      //const mic_m m_dist = lt(m_not_occluded,minDist,max_distance);
      sptr_node--;
      sptr_near--;
      
      if (unlikely(curNode == QBVH_TERMINAL_TOKEN)) break;
      STAT_COLLECTOR(StatCollector::chunk.numStackPops++);
      
      if (unlikely( m_dist == 0)) { continue;  }
      assert(curNode != QBVH_TERMINAL_TOKEN);
      
      while(1)
      {
        if (unlikely(qbvhLeaf(curNode))) break;
        assert (curNode != (unsigned int)(1 << 31));
        sptr_node--;
        sptr_near--;
        
        const Node* __restrict__ const qptr = qbvhChildPtr(qbvh,curNode);
        
        STAT_COLLECTOR(StatCollector::chunk.numTraversalSteps++);
        STAT_COLLECTOR(unsigned int quadhits=0);
        STAT_COLLECTOR(StatCollector::chunk.numActiveRaysTraversal+=countbits(m_dist));
        
        curNode = *sptr_node;
        minDist = *sptr_near;
        
        prefetch<PFHINT_L1>((mic_f*)qptr + 0); 
        prefetch<PFHINT_L1>((mic_f*)qptr + 1); 
        
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
          
          const mic_m lhit   = le(m_not_occluded,_max(min_distance,lnearP), _min(lfarP,max_distance));
          const mic_f boxDist = sel(lhit,lnearP,inf);
          
          const mic_m m_update_min = lt(m_not_occluded,boxDist,minDist);
          
          if (likely(nz_mask(lhit) != 0))
          {
            const unsigned int dd = qptr->lower[i].child;
            sptr_node++;
            sptr_near++;
            
            if (nz_mask(m_update_min) != 0) {
              *(sptr_node-1) = curNode;
              *(sptr_near-1) = minDist; 
              minDist = boxDist;
              curNode = dd;
            }
            else {
              *(sptr_node-1) = dd;
              *(sptr_near-1) = boxDist; 
            }
            assert(sptr_node - stack_node < QUAD_BVH_MAX_STACK_DEPTH);
          }	      
        }
        STAT_COLLECTOR(StatCollector::chunk.histo[quadhits]++);
      }
      
      assert (curNode != (unsigned int)(1 << 31));
      
      if (likely(curNode == QBVH_TERMINAL_TOKEN)) break;
      
      STAT_COLLECTOR(StatCollector::chunk.numActiveRaysAtLeaf+=countbits(lt(minDist,max_distance)));
      STAT_COLLECTOR(StatCollector::chunk.numLeafIntersections++);
      
      mic_m m_active_leaf = lt(m_not_occluded,minDist,max_distance);
      
      unsigned int itemOffset    = qbvhItemOffset(curNode);
      const Triangle1  *__restrict__ tptr    = (Triangle1*)((char*)accel + itemOffset);
      const unsigned int items = qbvhItems(curNode);
      
      assert(items > 0);
      assert(items <= 4);
      
      prefetch<PFHINT_L1>((mic_f*)tptr +  0); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  1); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  2); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  3); 
      
      STAT_COLLECTOR(StatCollector::chunk.numActiveRaysIntersection+=items * countbits(lt(minDist,max_distance)));
      
      for (unsigned int i=0;i<items;i++,tptr++) 
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
        
        //const mic_f normal = lcross_zxy(e1,e2);
        const mic_f normal = lcross_zxy(e2,e1);
        const mic_f nz = swAAAA(normal);
        const mic_f nx = swBBBB(normal);
        const mic_f ny = swCCCC(normal);
        const mic_f nom = madd231(madd231(v0_org_z * nz,v0_org_y, ny),v0_org_x,nx);
        const mic_f den = madd231(madd231(direction_z * nz,direction_y,ny),direction_x,nx); 
        
        const mic_f odz = msubr231(direction_x*v0_org_y,v0_org_x,direction_y);
        const mic_f odx = msubr231(direction_y*v0_org_z,v0_org_y,direction_z);
        const mic_f ody = msubr231(direction_z*v0_org_x,v0_org_z,direction_x);	  
        
        const mic_f u = madd231(madd231(odz*swCCCC(e2),ody,swBBBB(e2)),odx,swAAAA(e2));
        const mic_f v = madd231(madd231(odz*swCCCC(e1),ody,swBBBB(e1)),odx,swAAAA(e1));
        
        const mic_m den_gt_0 = gz(den); 
        
        const mic_m valid_u = (den_gt_0^ge(u,zero));
        const mic_m valid_v = (den_gt_0^ge(v,zero));
        const mic_m valid_uv = den_gt_0^ge(add_neg(u,v),den); 
        const mic_m m_aperture = m_not_occluded & valid_u & valid_v & valid_uv;
        
        const mic_f rcp_den = rcp(den);
        if ( unlikely(m_aperture == 0) ) continue;
        
        const mic_f t = rcp_den*nom;
        
#if defined(ENABLE_LIGHTSOURCES_MASK)
        const mic_i lsMask = mic_i(tptr->data[0]);
        mic_m m_final = lt(lt(m_aperture,min_distance,t),t,max_distance);
        m_final       = mask_test(m_final,lsMask,ray.data);
#else
        const mic_m m_final  = lt(lt(m_aperture,min_distance,t),t,max_distance);
#endif
        m_not_occluded &= ~m_final;
        m_active_leaf  &= ~m_final;
        if (unlikely(m_active_leaf == 0)) break;
      }
      if (unlikely(m_not_occluded == 0)) {
        return m_active;
      }
    }
    return m_active & ~m_not_occluded;
  }
  
  void BVH4AOSTriangle1Intersector16ChunkRefRegister () {
    TriangleMesh::intersectors16.add("bvh4aos","triangle1","chunk_ref","moeller",true ,BVH4AOSTriangle1Intersector16ChunkRef::create);
  }
}


