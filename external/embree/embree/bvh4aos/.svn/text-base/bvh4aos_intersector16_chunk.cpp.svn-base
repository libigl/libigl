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

#include "bvh4aos_intersector16_chunk.h"
#include "../geometry/triangles.h"

namespace embree
{
  template<typename TriangleIntersector>
  void BVH4AOSIntersector16Chunk<TriangleIntersector>::intersect(const BVH4AOSIntersector16Chunk* This, Ray16& ray, const __mmask valid_i)
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
      unsigned curNode = *(sptr_node-1);
      sptr_node--;
      sptr_near--;
      if (unlikely(curNode == BVH4AOS_INVALID_NODE)) 
        break;

      /* cull node if behind closest hit point */
      if (unlikely(none(m_dist))) continue;

      while (1)
      {
        /* test if this is a leaf node */
        if (unlikely(curNode & leaf_mask)) 
          break;
        
        /* pop of next node */
        sptr_node--;
        sptr_near--;
        const Node* const __restrict__ node = (Node*)((char*)nodes + (unsigned long)curNode);
        prefetch<PFHINT_L1>((mic_f*)node + 1); // depth first order prefetch		       
        curNode = *sptr_node;
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
              *(sptr_node-1) = curNode;
              *(sptr_near-1) = minDist; 
              minDist = childDist;
              curNode = child;
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
      if (unlikely(curNode == BVH4AOS_INVALID_NODE)) 
        break;
      
      /* decode leaf node */
      const unsigned itemOffset = BVH4AOS::itemOffset(curNode);
      const unsigned items      = BVH4AOS::itemCount (curNode);
      const Triangle* triangle = (Triangle*)((char*)triangles + itemOffset);
     
      /* intersect triangles */
      TriangleIntersector::intersect(valid, ray, triangle, items, bvh->vertices);
      max_distance = ray.tfar;
    }
  }

  template<typename TriangleIntersector>
  __mmask BVH4AOSIntersector16Chunk<TriangleIntersector>::occluded(const BVH4AOSIntersector16Chunk* This, Ray16& ray, const __mmask valid_i)
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
    
    mic_m m_not_occluded = valid;

    while (1)
    {
      /* pop next node from stack */
      mic_f minDist = upconv16f((const float*)(sptr_near-1));
      NodeRef curNode = *(sptr_node-1);
      sptr_node--;
      sptr_near--;
      if (unlikely(curNode.isInvalid())) 
        break;

      /* cull node if behind closest hit point */
      const mic_m m_dist = gt(m_not_occluded,max_distance,minDist);
      if (unlikely(none(m_dist))) continue;

      while (1)
      {
        /* test if this is a leaf node */
        if (unlikely(curNode.isLeaf())) 
          break;
        
        /* pop of next node */
        sptr_node--;
        sptr_near--;
        const Node* const node = curNode.node(nodes);
        prefetch<PFHINT_L1>((mic_f*)node + 1); // depth first order prefetch		
        curNode = *sptr_node;
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
              *(sptr_node-1) = curNode;
              *(sptr_near-1) = minDist; 
              minDist = childDist;
              curNode = child;
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
      if (unlikely(curNode.isInvalid())) 
        break;
      
       /* decode leaf node */
      size_t items;
      const Triangle* triangle = (Triangle*) curNode.leaf(triangles,items);
     
      /* intersect triangles */
      mic_m valid_leaf = lt(m_not_occluded,minDist,max_distance);
      mic_m mask = TriangleIntersector::occluded(valid_leaf, ray, triangle, items, bvh->vertices);

      /* update occlusion mask */
      m_not_occluded &= ~mask;
      if (unlikely(m_not_occluded == 0)) {
        return valid;
      }
    }
    return valid & ~m_not_occluded;
  }

  void BVH4AOSIntersector16ChunkRegister () 
  {
    //TriangleMesh::intersectors16.add("bvh4aos","triangle1i","chunk","moeller" ,false,BVH4AOSIntersector16Chunk<Triangle1iIntersector16<Intersector16MoellerTrumbore> >::create);
    //TriangleMesh::intersectors16.add("bvh4aos","triangle1i","chunk","pluecker",true ,BVH4AOSIntersector16Chunk<Triangle1iIntersector16<Intersector16Pluecker> >::create);
    //TriangleMesh::intersectors16.add("bvh4aos","triangle1v","chunk","moeller" ,false,BVH4AOSIntersector16Chunk<Triangle1vIntersector16<Intersector16MoellerTrumbore> >::create);
    //TriangleMesh::intersectors16.add("bvh4aos","triangle1v","chunk","pluecker",true ,BVH4AOSIntersector16Chunk<Triangle1vIntersector16<Intersector16Pluecker> >::create);
    //TriangleMesh::intersectors16.add("bvh4aos","triangle1v","chunk_ref","moeller",false ,BVH4AOSIntersector16Chunk<Triangle1vIntersector16MoellerTrumboreRefMIC>::create);
    //TriangleMesh::intersectors16.add("bvh4aos","triangle4i","chunk","moeller" ,false,BVH4AOSIntersector16Chunk<Triangle4iIntersector16<Intersector16MoellerTrumbore> >::create);
    //TriangleMesh::intersectors16.add("bvh4aos","triangle4i","chunk","pluecker",true ,BVH4AOSIntersector16Chunk<Triangle4iIntersector16<Intersector16Pluecker> >::create);
    //TriangleMesh::intersectors16.add("bvh4aos","triangle4v","chunk","moeller" ,false,BVH4AOSIntersector16Chunk<Triangle4vIntersector16<Intersector16MoellerTrumbore> >::create);
    //TriangleMesh::intersectors16.add("bvh4aos","triangle4v","chunk","pluecker",true ,BVH4AOSIntersector16Chunk<Triangle4vIntersector16<Intersector16Pluecker> >::create);
    //TriangleMesh::intersectors16.add("bvh4aos","triangle4aos","chunk","moeller",true ,BVH4AOSIntersector16Chunk<Triangle4AOSIntersector16MoellerTrumboreMIC>::create);
    TriangleMesh::intersectors16.add("bvh4aos","triangle1","chunk","moeller",true ,BVH4AOSIntersector16Chunk<Triangle1Intersector16MoellerTrumboreMIC>::create);
    //TriangleMesh::intersectors16.add("bvh4aos","triangle4" ,"chunk","moeller" ,true ,BVH4AOSIntersector16Chunk<Triangle4Intersector16MoellerTrumbore>::create);
    //TriangleMesh::intersectors16.add("bvh4aos","triangle8" ,"chunk","moeller" ,true ,BVH4AOSIntersector16Chunk<Triangle8Intersector16MoellerTrumbore>::create);

    VirtualScene::intersectors16.add("bvh4aos","virtual","chunk","virtual",true,BVH4AOSIntersector16Chunk<VirtualObjectIntersector16>::create);
    VirtualScene::intersectors16.setAccelDefaultTraverser("bvh4aos","chunk");
  }
}


