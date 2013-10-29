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

#include "bvh2_intersector4.h"
#include "geometry/triangles.h"

namespace embree
{
  /* ray/box intersection */
  __forceinline size_t intersectBox(const sse3f& org, const sse3f& rdir, const ssef& tnear, const ssef& tfar, const BVH2::Node* node, const int i, ssef& dist) 
  {
    const ssef dminx = (ssef(node->lower_upper_x[i+0]) - org.x) * rdir.x;
    const ssef dminy = (ssef(node->lower_upper_y[i+0]) - org.y) * rdir.y;
    const ssef dminz = (ssef(node->lower_upper_z[i+0]) - org.z) * rdir.z;
    const ssef dmaxx = (ssef(node->lower_upper_x[i+2]) - org.x) * rdir.x;
    const ssef dmaxy = (ssef(node->lower_upper_y[i+2]) - org.y) * rdir.y;
    const ssef dmaxz = (ssef(node->lower_upper_z[i+2]) - org.z) * rdir.z;
    
    const ssef dlowerx = min(dminx,dmaxx);
    const ssef dlowery = min(dminy,dmaxy);
    const ssef dlowerz = min(dminz,dmaxz);

    const ssef dupperx = max(dminx,dmaxx);
    const ssef duppery = max(dminy,dmaxy);
    const ssef dupperz = max(dminz,dmaxz);
    
    const ssef near = max(dlowerx,dlowery,dlowerz,tnear);
    const ssef far  = min(dupperx,duppery,dupperz,tfar );
    dist = near;
    
    return movemask(near <= far);
  }

  template<typename TriangleIntersector>
  void BVH2Intersector4Chunk<TriangleIntersector>::intersect(const BVH2Intersector4Chunk* This, Ray4& ray, const __m128 valid_i)
  {
    sseb valid = valid_i;
    const BVH2* bvh = This->bvh;
    STAT3(normal.travs,1,popcnt(valid),4);

    struct StackItem {
      NodeRef ptr;
      ssef dist;
    };
    
    StackItem stack[1+BVH2::maxDepth];   //!< stack of nodes that still need to get traversed
    StackItem* stackPtr = stack;                       //!< current stack pointer
    NodeRef cur = bvh->root;                        //!< in cur we track the ID of the current node

    /* let inactive rays miss all boxes */
    const sse3f rdir = rcp_safe(ray.dir);
    ray.tfar = select(valid,ray.tfar,ssef(neg_inf));

    while (true)
    {
      /*! downtraversal loop */
      while (likely(cur.isNode()))
      {
        STAT3(normal.trav_nodes,1,popcnt(valid),4);

        /* intersect packet with box of both children */
        const Node* node = cur.node();
        ssef dist0; size_t hit0 = intersectBox(ray.org,rdir,ray.tnear,ray.tfar,node,0,dist0);
        ssef dist1; size_t hit1 = intersectBox(ray.org,rdir,ray.tnear,ray.tfar,node,1,dist1);

        /*! if two children hit, push far node onto stack and continue with closer node */
        if (likely(hit0 != 0 && hit1 != 0)) {
          if (any(valid & (dist0 < dist1))) { 
            stackPtr->ptr = node->child(1); stackPtr->dist = dist1; stackPtr++; cur = node->child(0); 
          } else { 
            stackPtr->ptr = node->child(0); stackPtr->dist = dist0; stackPtr++; cur = node->child(1); 
          }
        }

        /*! if one child hit, continue with that child */
        else {
          if      (likely(hit0 != 0)) cur = node->child(0);
          else if (likely(hit1 != 0)) cur = node->child(1);
          else goto pop_node;
        }
      }

      /*! leaf node, intersect all triangles */
      {
        STAT3(normal.trav_leaves,1,popcnt(valid),4);
        size_t num; Triangle* tri = (Triangle*) cur.leaf(NULL,num);
        for (size_t i=0; i<num; i++) {
          TriangleIntersector::intersect(valid,ray,tri[i],bvh->vertices);
        }
      }

      /*! pop next node from stack */
pop_node:
      if (unlikely(stackPtr == stack)) break;
      --stackPtr;
      cur = stackPtr->ptr;
      if (unlikely(all(stackPtr->dist > ray.tfar))) goto pop_node;
    }
  }

  /* ray/box intersection */
  __forceinline size_t intersectBox(const sse3f& org, const sse3f& rdir, const ssef& tnear, const ssef& tfar, const BVH2::Node* node, const int i) 
  {
    const ssef dminx = (ssef(node->lower_upper_x[i+0]) - org.x) * rdir.x;
    const ssef dminy = (ssef(node->lower_upper_y[i+0]) - org.y) * rdir.y;
    const ssef dminz = (ssef(node->lower_upper_z[i+0]) - org.z) * rdir.z;
    const ssef dmaxx = (ssef(node->lower_upper_x[i+2]) - org.x) * rdir.x;
    const ssef dmaxy = (ssef(node->lower_upper_y[i+2]) - org.y) * rdir.y;
    const ssef dmaxz = (ssef(node->lower_upper_z[i+2]) - org.z) * rdir.z;
    
    const ssef dlowerx = min(dminx,dmaxx);
    const ssef dlowery = min(dminy,dmaxy);
    const ssef dlowerz = min(dminz,dmaxz);

    const ssef dupperx = max(dminx,dmaxx);
    const ssef duppery = max(dminy,dmaxy);
    const ssef dupperz = max(dminz,dmaxz);
    
    const ssef near = max(dlowerx,dlowery,dlowerz,tnear);
    const ssef far  = min(dupperx,duppery,dupperz,tfar );
    
    return movemask(near <= far);
  }

  template<typename TriangleIntersector>
  __m128 BVH2Intersector4Chunk<TriangleIntersector>::occluded(const BVH2Intersector4Chunk* This, Ray4& ray, const __m128 valid_i)
  {
    sseb valid = valid_i;
    sseb terminated = !valid;
    const BVH2* bvh = This->bvh;
    STAT3(shadow.travs,1,popcnt(valid),4);

    NodeRef stack[1+BVH2::maxDepth]; //!< stack of nodes that still need to get traversed
    NodeRef* stackPtr = stack;                    //!< current stack pointer
    NodeRef cur = bvh->root;                      //!< in cur we track the ID of the current node

    /* let inactive rays miss all boxes */
    const sse3f rdir = rcp_safe(ray.dir);
    ssef rayFar  = select(terminated,ssef(neg_inf),ray.tfar);

    while (true)
    {
      /*! downtraversal loop */
      while (likely(cur.isNode()))
      {
        STAT3(normal.trav_nodes,1,popcnt(valid),4);

        /* intersect packet with box of both children */
        const Node* node = cur.node();
        const size_t hit0 = intersectBox(ray.org,rdir,ray.tnear,rayFar,node,0);
        const size_t hit1 = intersectBox(ray.org,rdir,ray.tnear,rayFar,node,1);
        
        /*! if two children are hit push both onto stack */
        if (likely(hit0 != 0 && hit1 != 0)) {
          *stackPtr = node->child(0); stackPtr++; cur = node->child(1);
        }
        
        /*! if one child hit, continue with that child */
        else {
          if      (likely(hit0 != 0)) cur = node->child(0);
          else if (likely(hit1 != 0)) cur = node->child(1);
          else goto pop_node;
        }
      }

      /*! leaf node, intersect all triangles */
      {
        STAT3(shadow.trav_leaves,1,popcnt(valid),4);
        size_t num; Triangle* tri = (Triangle*) cur.leaf(NULL,num);
        for (size_t i=0; i<num; i++) {
          terminated |= TriangleIntersector::occluded(valid,ray,tri[i],bvh->vertices);
          if (all(terminated)) return terminated;
        }

        /* let terminated rays miss all boxes */
        rayFar = select(terminated,ssef(neg_inf),rayFar);
      }

      /*! pop next node from stack */
pop_node:
      if (unlikely(stackPtr == stack)) break;
      cur = *(--stackPtr);
    }
    return terminated;
  }

  void BVH2Intersector4ChunkRegister () 
  {
    TriangleMesh::intersectors4.add("bvh2","triangle1i","chunk","moeller" ,false,BVH2Intersector4Chunk<Triangle1iIntersector4<Intersector4MoellerTrumbore> >::create);
    TriangleMesh::intersectors4.add("bvh2","triangle1i","chunk","pluecker",true ,BVH2Intersector4Chunk<Triangle1iIntersector4<Intersector4Pluecker> >::create);
    TriangleMesh::intersectors4.add("bvh2","triangle1v","chunk","moeller" ,false,BVH2Intersector4Chunk<Triangle1vIntersector4<Intersector4MoellerTrumbore> >::create);
    TriangleMesh::intersectors4.add("bvh2","triangle1v","chunk","pluecker",true ,BVH2Intersector4Chunk<Triangle1vIntersector4<Intersector4Pluecker> >::create);
    TriangleMesh::intersectors4.add("bvh2","triangle4i","chunk","moeller" ,false,BVH2Intersector4Chunk<Triangle4iIntersector4<Intersector4MoellerTrumbore> >::create);
    TriangleMesh::intersectors4.add("bvh2","triangle4i","chunk","pluecker",true ,BVH2Intersector4Chunk<Triangle4iIntersector4<Intersector4Pluecker> >::create);
    TriangleMesh::intersectors4.add("bvh2","triangle4v","chunk","moeller" ,false,BVH2Intersector4Chunk<Triangle4vIntersector4<Intersector4MoellerTrumbore> >::create);
    TriangleMesh::intersectors4.add("bvh2","triangle4v","chunk","pluecker",true ,BVH2Intersector4Chunk<Triangle4vIntersector4<Intersector4Pluecker> >::create);
    TriangleMesh::intersectors4.add("bvh2","triangle1" ,"chunk","moeller" ,true ,BVH2Intersector4Chunk<Triangle1Intersector4MoellerTrumbore>::create);
    TriangleMesh::intersectors4.add("bvh2","triangle4" ,"chunk","moeller" ,true ,BVH2Intersector4Chunk<Triangle4Intersector4MoellerTrumbore>::create);
    //TriangleMesh::intersectors4.add("bvh2","triangle8" ,"chunk","moeller" ,true ,BVH2Intersector4Chunk<Triangle8Intersector4MoellerTrumbore>::create);
    TriangleMesh::intersectors4.setAccelDefaultTraverser("bvh2","chunk");

    VirtualScene::intersectors4.add("bvh2","virtual"   ,"chunk","virtual" ,true ,BVH2Intersector4Chunk<VirtualObjectIntersector4>::create);
    VirtualScene::intersectors4.setAccelDefaultTraverser("bvh2","chunk");
  }
}
