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

#include "bvh4_intersector4_chunk.h"
#include "../geometry/triangles.h"

namespace embree
{
  template<typename TriangleIntersector4>
  void BVH4Intersector4Chunk<TriangleIntersector4>::intersect(const BVH4Intersector4Chunk* This, Ray4& ray, const __m128 valid_i)
  {
    sseb valid = valid_i;
    NodeRef invalid = (NodeRef)1;
    const BVH4* bvh = This->bvh;
    STAT3(normal.travs,1,popcnt(valid),4);

    /* load ray into registers */
    ssef ray_near = select(valid,ray.tnear,pos_inf);
    ssef ray_far  = select(valid,ray.tfar ,neg_inf);
    sse3f rdir = rcp_safe(ray.dir);
    ray.tfar = ray_far;

    /* allocate stack and push root node */
    NodeRef stack_node[3*BVH4::maxDepth+1];
    ssef  stack_near[3*BVH4::maxDepth+1];
    stack_node[0] = invalid;
    stack_near[0] = ssef(inf);
    stack_node[1] = bvh->root;
    stack_near[1] = ray_near;
    NodeRef* sptr_node = stack_node+2;
    ssef * sptr_near = stack_near+2;
 
    while (1)
    {
      /* pop next node from stack */
      sptr_node--;
      sptr_near--;
      ssef  curDist = *sptr_near;
      NodeRef curNode = *sptr_node;
      if (unlikely(curNode == invalid))
        break;

      /* cull node if behind closest hit point */
      const sseb m_dist = curDist < ray_far;
      if (unlikely(none(m_dist))) 
        continue;

      while (1)
      {
        /* test if this is a leaf node */
        if (unlikely(curNode.isLeaf())) 
          break;

        STAT3(normal.trav_nodes,1,popcnt(valid),4);
        
        const Node* const node = curNode.node(bvh->nodePtr()); //NodeRef(curNode).node(nodes);
        //prefetch<PFHINT_L1>((ssef*)node + 1); // depth first order prefetch	
        
        /* pop of next node */
        sptr_node--;
        sptr_near--;
        curNode = *sptr_node;
        curDist = *sptr_near;
                
        for (unsigned i=0;i<4;i++)
	{
          const ssef dminx = (ssef(node->lower_x[i]) - ray.org.x) * rdir.x;
          const ssef dmaxx = (ssef(node->upper_x[i]) - ray.org.x) * rdir.x;
          const ssef dminy = (ssef(node->lower_y[i]) - ray.org.y) * rdir.y;
          const ssef dmaxy = (ssef(node->upper_y[i]) - ray.org.y) * rdir.y;
          const ssef dminz = (ssef(node->lower_z[i]) - ray.org.z) * rdir.z;
          const ssef dmaxz = (ssef(node->upper_z[i]) - ray.org.z) * rdir.z;
          
          const ssef dlowerx = min(dminx,dmaxx);
          const ssef dupperx = max(dminx,dmaxx);
          const ssef dlowery = min(dminy,dmaxy);
          const ssef duppery = max(dminy,dmaxy);
          const ssef dlowerz = min(dminz,dmaxz);
          const ssef dupperz = max(dminz,dmaxz);
          
          const ssef near = max(dlowerx,dlowery,dlowerz,ray_near);
          const ssef far  = min(dupperx,duppery,dupperz,ray_far );
          const sseb mhit = near <= far;
          
          const ssef childDist = select(mhit,near,inf);
          const sseb closer = childDist < curDist;

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
      STAT3(normal.trav_leaves,1,popcnt(valid),4);
      Triangle* tri = (Triangle*) curNode.leaf(bvh->triPtr(),num);
      
      /* intersect triangles */
      for (size_t i=0; i<num; i++)
        TriangleIntersector4::intersect(valid,ray,tri[i],bvh->vertices);
      
      ray_far = ray.tfar;
    }
  }
    
  template<typename TriangleIntersector4>
  __m128 BVH4Intersector4Chunk<TriangleIntersector4>::occluded(const BVH4Intersector4Chunk* This, Ray4& ray, const __m128 valid_i) 
  {
    sseb valid = valid_i;
    NodeRef invalid = (NodeRef)1;
    sseb terminated = !valid;
    const BVH4* bvh = This->bvh;
    STAT3(shadow.travs,1,popcnt(valid),4);

    /* load ray into registers */
    ssef ray_near = select(valid,ray.tnear,pos_inf);
    ssef ray_far  = select(valid,ray.tfar ,neg_inf);
    sse3f rdir = rcp_safe(ray.dir);

    /* allocate stack and push root node */
    NodeRef stack_node[3*BVH4::maxDepth+1];
    ssef  stack_near[3*BVH4::maxDepth+1];
    stack_node[0] = invalid;
    stack_near[0] = ssef(inf);
    stack_node[1] = bvh->root;
    stack_near[1] = ray_near;
    NodeRef* sptr_node = stack_node+2;
    ssef * sptr_near = stack_near+2;
 
    while (1)
    {
      /* pop next node from stack */
      sptr_node--;
      sptr_near--;
      ssef  curDist = *sptr_near;
      NodeRef curNode = *sptr_node;
      if (unlikely(curNode == invalid))
        break;

      /* cull node if behind closest hit point */
      const sseb m_dist = curDist < ray_far;
      if (unlikely(none(m_dist))) 
        continue;

      while (1)
      {
        /* test if this is a leaf node */
        if (unlikely(curNode.isLeaf())) 
          break;
        
        STAT3(shadow.trav_nodes,1,popcnt(valid),4);
        const Node* const node = curNode.node(bvh->nodePtr()); //NodeRef(curNode).node(nodes);
        //prefetch<PFHINT_L1>((ssef*)node + 1); // depth first order prefetch	
        
        /* pop of next node */
        sptr_node--;
        sptr_near--;
        curNode = *sptr_node;
        curDist = *sptr_near;
                
        for (unsigned i=0;i<4;i++)
	{
          const ssef dminx = (ssef(node->lower_x[i]) - ray.org.x) * rdir.x;
          const ssef dmaxx = (ssef(node->upper_x[i]) - ray.org.x) * rdir.x;
          const ssef dminy = (ssef(node->lower_y[i]) - ray.org.y) * rdir.y;
          const ssef dmaxy = (ssef(node->upper_y[i]) - ray.org.y) * rdir.y;
          const ssef dminz = (ssef(node->lower_z[i]) - ray.org.z) * rdir.z;
          const ssef dmaxz = (ssef(node->upper_z[i]) - ray.org.z) * rdir.z;
          
          const ssef dlowerx = min(dminx,dmaxx);
          const ssef dupperx = max(dminx,dmaxx);
          const ssef dlowery = min(dminy,dmaxy);
          const ssef duppery = max(dminy,dmaxy);
          const ssef dlowerz = min(dminz,dmaxz);
          const ssef dupperz = max(dminz,dmaxz);
          
          const ssef near = max(dlowerx,dlowery,dlowerz,ray_near);
          const ssef far  = min(dupperx,duppery,dupperz,ray_far );
          const sseb mhit = near <= far;
          
          const ssef childDist = select(mhit,near,inf);
          const sseb closer = childDist < curDist;

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
      STAT3(shadow.trav_leaves,1,popcnt(valid),4);
      Triangle* tri = (Triangle*) curNode.leaf(bvh->triPtr(),num);
      
      /* intersect triangles */
      for (size_t i=0; i<num; i++) {
        terminated |= TriangleIntersector4::occluded(valid,ray,tri[i],bvh->vertices);
        if (all(terminated)) return terminated;
      }
      ray_far = select(terminated,neg_inf,ray_far);
    }
    return valid & terminated;
  }

  void BVH4Intersector4ChunkRegister () 
  {
    TriangleMesh::intersectors4.add("bvh4","triangle1i","chunk","moeller" ,false,BVH4Intersector4Chunk<Triangle1iIntersector4<Intersector4MoellerTrumbore> >::create);
    TriangleMesh::intersectors4.add("bvh4","triangle1i","chunk","pluecker",true ,BVH4Intersector4Chunk<Triangle1iIntersector4<Intersector4Pluecker> >::create);
    TriangleMesh::intersectors4.add("bvh4","triangle1v","chunk","moeller" ,false,BVH4Intersector4Chunk<Triangle1vIntersector4<Intersector4MoellerTrumbore> >::create);
    TriangleMesh::intersectors4.add("bvh4","triangle1v","chunk","pluecker",true ,BVH4Intersector4Chunk<Triangle1vIntersector4<Intersector4Pluecker> >::create);
    TriangleMesh::intersectors4.add("bvh4","triangle4i","chunk","moeller" ,false,BVH4Intersector4Chunk<Triangle4iIntersector4<Intersector4MoellerTrumbore> >::create);
    TriangleMesh::intersectors4.add("bvh4","triangle4i","chunk","pluecker",true ,BVH4Intersector4Chunk<Triangle4iIntersector4<Intersector4Pluecker> >::create);
    TriangleMesh::intersectors4.add("bvh4","triangle4v","chunk","moeller" ,false,BVH4Intersector4Chunk<Triangle4vIntersector4<Intersector4MoellerTrumbore> >::create);
    TriangleMesh::intersectors4.add("bvh4","triangle4v","chunk","pluecker",true ,BVH4Intersector4Chunk<Triangle4vIntersector4<Intersector4Pluecker> >::create);
    TriangleMesh::intersectors4.add("bvh4","triangle1" ,"chunk","moeller" ,true ,BVH4Intersector4Chunk<Triangle1Intersector4MoellerTrumbore>::create);
    TriangleMesh::intersectors4.add("bvh4","triangle4" ,"chunk","moeller" ,true ,BVH4Intersector4Chunk<Triangle4Intersector4MoellerTrumbore>::create);
    //TriangleMesh::intersectors4.add("bvh4","triangle8" ,"chunk","moeller" ,true ,BVH4Intersector4Chunk<Triangle8Intersector4MoellerTrumbore>::create);
    TriangleMesh::intersectors4.setAccelDefaultTraverser("bvh4","chunk");

    VirtualScene::intersectors4.add("bvh4","virtual"   ,"chunk","virtual" ,true ,BVH4Intersector4Chunk<VirtualObjectIntersector4>::create);
    VirtualScene::intersectors4.setAccelDefaultTraverser("bvh4","chunk");
  }
}
