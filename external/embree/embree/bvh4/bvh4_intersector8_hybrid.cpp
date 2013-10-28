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

#if __AVX__

#include "bvh4_intersector8_hybrid.h"
#include "bvh4_intersector1.h"
#include "../geometry/triangles.h"

#define SWITCH_THRESHOLD 3

namespace embree
{
  template<typename TriangleIntersector8, typename TriangleIntersector1>
  void BVH4Intersector8Hybrid<TriangleIntersector8, TriangleIntersector1>::intersect(const BVH4Intersector8Hybrid* This, Ray8& ray, const __m256 valid_i)
  {
    avxb valid = valid_i;
    const BVH4* bvh = This->bvh;
    STAT3(normal.travs,1,popcnt(valid),8);

    NodeRef invalid = (NodeRef)1;

    /* load ray into registers */
    avxf ray_near = select(valid,ray.tnear,pos_inf);
    avxf ray_far  = select(valid,ray.tfar ,neg_inf);
    avx3f rdir = rcp_safe(ray.dir);
    ray.tfar = ray_far;

    /* allocate stack and push root node */
    NodeRef stack_node[3*BVH4::maxDepth+1];
    avxf  stack_near[3*BVH4::maxDepth+1];
    stack_node[0] = invalid;
    stack_near[0] = inf;
    stack_node[1] = bvh->root;
    stack_near[1] = ray_near;
    NodeRef* sptr_node = stack_node+2;
    avxf * sptr_near = stack_near+2;
 
    while (1)
    {
      /* pop next node from stack */
      sptr_node--;
      sptr_near--;
      avxf  curDist = *sptr_near;
      NodeRef curNode = *sptr_node;
      if (unlikely(curNode == invalid))
        break;

      /* cull node if behind closest hit point */
      const avxb active = curDist < ray_far;
      if (unlikely(none(active))) 
        continue;

      /* switch to single ray traversal */
      size_t bits = movemask(active);
      if (unlikely(__popcnt(bits) <= SWITCH_THRESHOLD)) {
        for (size_t i=__bsf(bits); bits!=0; bits=__btc(bits,i), i=__bsf(bits)) {
          BVH4Intersector1<TriangleIntersector1>::intersect1(bvh,curNode,i,ray,rdir);
        }
        ray_far = ray.tfar;
        continue;
      }

      while (1)
      {
        /* test if this is a leaf node */
        if (unlikely(curNode.isLeaf())) 
          break;
        
        const Node* const node = curNode.node(bvh->nodePtr()); //NodeRef(curNode).node(nodes);
        
        /* pop of next node */
        sptr_node--;
        sptr_near--;
        curNode = *sptr_node;
        curDist = *sptr_near;
                
        for (unsigned i=0; i<4; i++)
	{
          const avxf dminx = (avxf(node->lower_x[i]) - ray.org.x) * rdir.x;
          const avxf dmaxx = (avxf(node->upper_x[i]) - ray.org.x) * rdir.x;
          const avxf dminy = (avxf(node->lower_y[i]) - ray.org.y) * rdir.y;
          const avxf dmaxy = (avxf(node->upper_y[i]) - ray.org.y) * rdir.y;
          const avxf dminz = (avxf(node->lower_z[i]) - ray.org.z) * rdir.z;
          const avxf dmaxz = (avxf(node->upper_z[i]) - ray.org.z) * rdir.z;
          const NodeRef child = node->child(i);
          
          const avxf dlowerx = min(dminx,dmaxx);
          const avxf dupperx = max(dminx,dmaxx);
          const avxf dlowery = min(dminy,dmaxy);
          const avxf duppery = max(dminy,dmaxy);
          const avxf dlowerz = min(dminz,dmaxz);
          const avxf dupperz = max(dminz,dmaxz);
          
          const avxf near = max(dlowerx,dlowery,dlowerz,ray_near);
          const avxf far  = min(dupperx,duppery,dupperz,ray_far );
          const avxb mhit = near <= far;
          
          const avxf childDist = select(mhit,near,inf);
          const avxb closer = childDist < curDist;

          /* if we hit the child we choose to continue with that child if it 
             is closer than the current next child, or we push it onto the stack */
          if (likely(any(mhit)))
          {
            //if (child != invalid)
            {
              sptr_node++;
              sptr_near++;
            
              /* push cur node onto stack and continue with hit child */
              if (any(closer)) {
                *(sptr_node-1) = curNode;
                *(sptr_near-1) = curDist; 
                curDist = childDist;
                curNode = child;
              } 
              
              /* push hit child onto stack*/
              else {
                *(sptr_node-1) = child;
                *(sptr_near-1) = childDist; 
              }          
            }	      
          }
        }
      }
      
      /* return if stack is empty */
      if (unlikely(curNode == invalid)) 
        break;
      
      /* decode leaf node */
      size_t num;
      Triangle* tri = (Triangle*) curNode.leaf(bvh->triPtr(),num);
      
      /* intersect triangles */
      for (size_t i=0; i<num; i++)
        TriangleIntersector8::intersect(valid,ray,tri[i],bvh->vertices);
      
      ray_far = ray.tfar;
    }
  }

  template<typename TriangleIntersector8, typename TriangleIntersector1>
  __m256 BVH4Intersector8Hybrid<TriangleIntersector8, TriangleIntersector1>::occluded(const BVH4Intersector8Hybrid* This, Ray8& ray, const __m256 valid_i)
  {
    avxb valid = valid_i;
    NodeRef invalid = (NodeRef)1;
    avxb terminated = !valid;
    const BVH4* bvh = This->bvh;
    
    /* load ray into registers */
    avxf ray_near = select(valid,ray.tnear,pos_inf);
    avxf ray_far  = select(valid,ray.tfar ,neg_inf);
    avx3f rdir = rcp_safe(ray.dir);

    /* allocate stack and push root node */
    NodeRef stack_node[3*BVH4::maxDepth+1];
    avxf  stack_near[3*BVH4::maxDepth+1];
    stack_node[0] = invalid;
    stack_near[0] = inf;
    stack_node[1] = bvh->root;
    stack_near[1] = ray_near;
    NodeRef* sptr_node = stack_node+2;
    avxf * sptr_near = stack_near+2;
 
    while (1)
    {
      /* pop next node from stack */
      sptr_node--;
      sptr_near--;
      avxf  curDist = *sptr_near;
      NodeRef curNode = *sptr_node;
      if (unlikely(curNode == invalid))
        break;

      /* cull node if behind closest hit point */
      const avxb active = curDist < ray_far;
      if (unlikely(none(active))) 
        continue;

      /* switch to single ray traversal */
      size_t bits = movemask(active);
      if (unlikely(__popcnt(bits) <= SWITCH_THRESHOLD)) {
        for (size_t i=__bsf(bits); bits!=0; bits=__btc(bits,i), i=__bsf(bits)) {
          if (BVH4Intersector1<TriangleIntersector1>::occluded1(bvh,curNode,i,ray,rdir))
            terminated[i] = -1;
        }
        if (all(terminated)) return terminated;

        /* let terminated rays miss all boxes */
        ray_far = select(terminated,neg_inf,ray_far);
        continue;
      }

      while (1)
      {
        /* test if this is a leaf node */
        if (unlikely(curNode.isLeaf())) 
          break;
        
        const Node* const node = curNode.node(bvh->nodePtr()); //NodeRef(curNode).node(nodes);
        //prefetch<PFHINT_L1>((avxf*)node + 1); // depth first order prefetch	
        
        /* pop of next node */
        sptr_node--;
        sptr_near--;
        curNode = *sptr_node;
        curDist = *sptr_near;
                
        for (unsigned i=0;i<4;i++)
	{
          const avxf dminx = (avxf(node->lower_x[i]) - ray.org.x) * rdir.x;
          const avxf dmaxx = (avxf(node->upper_x[i]) - ray.org.x) * rdir.x;
          const avxf dminy = (avxf(node->lower_y[i]) - ray.org.y) * rdir.y;
          const avxf dmaxy = (avxf(node->upper_y[i]) - ray.org.y) * rdir.y;
          const avxf dminz = (avxf(node->lower_z[i]) - ray.org.z) * rdir.z;
          const avxf dmaxz = (avxf(node->upper_z[i]) - ray.org.z) * rdir.z;
          
          const avxf dlowerx = min(dminx,dmaxx);
          const avxf dupperx = max(dminx,dmaxx);
          const avxf dlowery = min(dminy,dmaxy);
          const avxf duppery = max(dminy,dmaxy);
          const avxf dlowerz = min(dminz,dmaxz);
          const avxf dupperz = max(dminz,dmaxz);
          
          const avxf near = max(dlowerx,dlowery,dlowerz,ray_near);
          const avxf far  = min(dupperx,duppery,dupperz,ray_far );
          const avxb mhit = near <= far;
          
          const avxf childDist = select(mhit,near,inf);
          const avxb closer = childDist < curDist;

          /* if we hit the child we choose to continue with that child if it 
             is closer than the current next child, or we push it onto the stack */
          if (likely(any(mhit)))
          {
            const NodeRef child = node->child(i);
            //if (child != invalid)
            {
              sptr_node++;
              sptr_near++;
            
              /* push cur node onto stack and continue with hit child */
              if (any(closer)) {
                *(sptr_node-1) = curNode;
                *(sptr_near-1) = curDist; 
                curDist = childDist;
                curNode = child;
              } 
              
              /* push hit child onto stack*/
              else {
                *(sptr_node-1) = child;
                *(sptr_near-1) = childDist; 
              }          
            }	      
          }
        }
      }
      
      /* return if stack is empty */
      if (unlikely(curNode == invalid)) 
        break;
      
      /* decode leaf node */
      size_t num;
      Triangle* tri = (Triangle*) curNode.leaf(bvh->triPtr(),num);
      
      /* intersect triangles */
      for (size_t i=0; i<num; i++) {
        terminated |= TriangleIntersector8::occluded(valid,ray,tri[i],bvh->vertices);
        if (all(terminated)) return terminated;
      }
      ray_far = select(terminated,neg_inf,ray_far);
    }
    return valid & terminated;
  }

  void BVH4Intersector8HybridRegister () 
  {
    TriangleMesh::intersectors8.add("bvh4","triangle1i","hybrid","moeller" ,false, BVH4Intersector8Hybrid<Triangle1iIntersector8<Intersector8MoellerTrumbore>, Triangle1iIntersector1<Intersector1MoellerTrumbore> >::create);
    //TriangleMesh::intersectors8.add("bvh4","triangle1i","hybrid","pluecker",true ,BVH4Intersector8Hybrid<Triangle1iIntersector8<Intersector8Pluecker> >::create);
    //TriangleMesh::intersectors8.add("bvh4","triangle1v","hybrid","moeller" ,false,BVH4Intersector8Hybrid<Triangle1vIntersector8<Intersector8MoellerTrumbore> >::create);
    //TriangleMesh::intersectors8.add("bvh4","triangle1v","hybrid","pluecker",true ,BVH4Intersector8Hybrid<Triangle1vIntersector8<Intersector8Pluecker> >::create);
    //TriangleMesh::intersectors8.add("bvh4","triangle4i","hybrid","moeller" ,false,BVH4Intersector8Hybrid<Triangle4iIntersector8<Intersector8MoellerTrumbore> >::create);
    //TriangleMesh::intersectors8.add("bvh4","triangle4i","hybrid","pluecker",true ,BVH4Intersector8Hybrid<Triangle4iIntersector8<Intersector8Pluecker> >::create);
    //TriangleMesh::intersectors8.add("bvh4","triangle4v","hybrid","moeller" ,false,BVH4Intersector8Hybrid<Triangle4vIntersector8<Intersector8MoellerTrumbore> >::create);
    //TriangleMesh::intersectors8.add("bvh4","triangle4v","hybrid","pluecker",true ,BVH4Intersector8Hybrid<Triangle4vIntersector8<Intersector8Pluecker> >::create);
    // TriangleMesh::intersectors8.add("bvh4","triangle1" ,"hybrid","moeller" ,true ,BVH4Intersector8Hybrid<Triangle1Intersector8MoellerTrumbore>::create);
    TriangleMesh::intersectors8.add("bvh4","triangle4" ,"hybrid","moeller" ,true ,BVH4Intersector8Hybrid<Triangle4Intersector8MoellerTrumbore, Triangle4Intersector1MoellerTrumbore>::create);
    //TriangleMesh::intersectors8.add("bvh4","triangle8" ,"hybrid","moeller" ,true ,BVH4Intersector8Hybrid<Triangle8Intersector8MoellerTrumbore>::create);
    TriangleMesh::intersectors8.setAccelDefaultTraverser("bvh4.triangle4","hybrid");

    //VirtualScene::intersectors8.add("bvh4","virtual","hybrid","virtual" ,true,BVH4Intersector8Hybrid<VirtualObjectIntersector8>::create);
    //VirtualScene::intersectors8.setAccelDefaultTraverser("bvh4","hybrid");
  }
}
#endif

