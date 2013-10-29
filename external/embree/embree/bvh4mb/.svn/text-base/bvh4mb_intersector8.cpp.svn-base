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

#include "bvh4mb_intersector8.h"
#include "../geometry/triangles.h"

namespace embree
{
  /*! An item on the stack holds the node and distance. */
  struct StackItemBVH4MBPacket8 
  {
    /*! pointer to the node */
    void* ptr;

    /*! distance of stack node */
    float dist;

    /*! Sort 3 stack items. */
    __forceinline friend void sort(StackItemBVH4MBPacket8& s1, StackItemBVH4MBPacket8& s2, StackItemBVH4MBPacket8& s3)
    {
      if (s2.dist < s1.dist) std::swap(s2,s1);
      if (s3.dist < s2.dist) std::swap(s3,s2);
      if (s2.dist < s1.dist) std::swap(s2,s1);
    }
    
    /*! Sort 4 stack items. */
    __forceinline friend void sort(StackItemBVH4MBPacket8& s1, StackItemBVH4MBPacket8& s2, StackItemBVH4MBPacket8& s3, StackItemBVH4MBPacket8& s4)
    {
      if (s2.dist < s1.dist) std::swap(s2,s1);
      if (s4.dist < s3.dist) std::swap(s4,s3);
      if (s3.dist < s1.dist) std::swap(s3,s1);
      if (s4.dist < s2.dist) std::swap(s4,s2);
      if (s3.dist < s2.dist) std::swap(s3,s2);
    }
  };

  /* ray/box intersection */
  __forceinline size_t intersectBox(const Ray8& ray, const avx3f& rdir, const BVH4MB::Node* node, const int i, avxf& dist) 
  {
    if (unlikely(node->child[i]->isEmptyLeaf())) return 0;

    const avxf lower_x = avxf(node->lower_x[i]) + ray.time * avxf(node->lower_dx[i]);
    const avxf lower_y = avxf(node->lower_y[i]) + ray.time * avxf(node->lower_dy[i]);
    const avxf lower_z = avxf(node->lower_z[i]) + ray.time * avxf(node->lower_dz[i]);
    const avxf upper_x = avxf(node->upper_x[i]) + ray.time * avxf(node->upper_dx[i]);
    const avxf upper_y = avxf(node->upper_y[i]) + ray.time * avxf(node->upper_dy[i]);
    const avxf upper_z = avxf(node->upper_z[i]) + ray.time * avxf(node->upper_dz[i]);

    const avxf dminx = (lower_x - ray.org.x) * rdir.x;
    const avxf dminy = (lower_y - ray.org.y) * rdir.y;
    const avxf dminz = (lower_z - ray.org.z) * rdir.z;
    const avxf dmaxx = (upper_x - ray.org.x) * rdir.x;
    const avxf dmaxy = (upper_y - ray.org.y) * rdir.y;
    const avxf dmaxz = (upper_z - ray.org.z) * rdir.z;
    
    const avxf dlowerx = min(dminx,dmaxx);
    const avxf dlowery = min(dminy,dmaxy);
    const avxf dlowerz = min(dminz,dmaxz);
    
    const avxf dupperx = max(dminx,dmaxx);
    const avxf duppery = max(dminy,dmaxy);
    const avxf dupperz = max(dminz,dmaxz);
    
    const avxf near = max(dlowerx,dlowery,dlowerz,ray.tnear);
    const avxf far  = min(dupperx,duppery,dupperz,ray.tfar );
    dist = near;
    
    return movemask(near <= far);
  }

  template<typename TriangleIntersector>
  void BVH4MBIntersector8Chunk<TriangleIntersector>::intersect(const BVH4MBIntersector8Chunk* This, Ray8& ray, const __m256 valid_i)
  {
    avxb valid = valid_i;
    const BVH4MB* bvh = This->bvh;
    STAT3(normal.travs,1,popcnt(valid),4);
    
    StackItemBVH4MBPacket8 stack[2+3*BVH4MB::maxDepth];
    StackItemBVH4MBPacket8* stackPtr = stack+1; //!< current stack pointer
    stack[0].ptr = bvh->root; 
    stack[0].dist = neg_inf;

    /* let inactive rays miss all boxes */
    avx3f rdir = rcp_safe(ray.dir);
    ray.tfar = select(valid,ray.tfar,avxf(neg_inf));

    while (true)
    {
      /*! pop next node */
      if (unlikely(stackPtr == stack)) break;
      stackPtr--;
      BVH4MB::Base* cur  = (BVH4MB::Base*) stackPtr->ptr;
 
      /* if popped node if too far, pop next one */
      if (unlikely(all(stackPtr->dist > ray.tfar)))
        continue;

    next:
        
      /* this is an inner node */
      if (likely(cur->isNode()))
      {
        STAT3(normal.trav_nodes,1,popcnt(valid),4);

        /* intersect packet with all boxes */
        const BVH4MB::Node* node = cur->node();
        avxf dist0; size_t hit0 = intersectBox(ray,rdir,node,0,dist0);
        avxf dist1; size_t hit1 = intersectBox(ray,rdir,node,1,dist1);
        avxf dist2; size_t hit2 = intersectBox(ray,rdir,node,2,dist2);
        avxf dist3; size_t hit3 = intersectBox(ray,rdir,node,3,dist3);
        
        /* push hit nodes onto stack */
        size_t cnt = 0;
        if (unlikely(hit0 != 0)) { stackPtr->ptr = node->child[0]; stackPtr->dist = dist0[__bsf(hit0)]; stackPtr++; cnt++; }
        if (unlikely(hit1 != 0)) { stackPtr->ptr = node->child[1]; stackPtr->dist = dist1[__bsf(hit1)]; stackPtr++; cnt++; }
        if (unlikely(hit2 != 0)) { stackPtr->ptr = node->child[2]; stackPtr->dist = dist2[__bsf(hit2)]; stackPtr++; cnt++; }
        if (unlikely(hit3 != 0)) { stackPtr->ptr = node->child[3]; stackPtr->dist = dist3[__bsf(hit3)]; stackPtr++; cnt++; }
        
        /*! if no child is hit, pop next node */
        if (unlikely(cnt == 0))
          continue;

        else if (unlikely(cnt == 1)) {
          cur = (BVH4MB::Base*) stackPtr[-1].ptr; stackPtr--;
          goto next;
        }
        
        else if (unlikely(cnt == 2)) {
          if (likely(stackPtr[-2].dist >= stackPtr[-1].dist)) {
            cur = (BVH4MB::Base*) stackPtr[-1].ptr; stackPtr--;
            goto next;
          }
          std::swap(stackPtr[-2],stackPtr[-1]);
          cur = (BVH4MB::Base*) stackPtr[-1].ptr; stackPtr--;
          goto next;
        }
        else if (unlikely(cnt == 3)) {
          sort(stackPtr[-1],stackPtr[-2],stackPtr[-3]);
          cur = (BVH4MB::Base*) stackPtr[-1].ptr; stackPtr--;
          goto next;
        } 
        else if (unlikely(cnt == 4)) {
          sort(stackPtr[-1],stackPtr[-2],stackPtr[-3],stackPtr[-4]);
          cur = (BVH4MB::Base*) stackPtr[-1].ptr; stackPtr--;
          goto next;
        } 
      } 

      /* this is a leaf node */
      else 
      {
        STAT3(normal.trav_leaves,1,popcnt(valid),4);
        size_t num; Triangle* tri = (Triangle*) cur->leaf(num);
        for (size_t i=0; i<num; i++)
          TriangleIntersector::intersect(valid,ray,tri[i],bvh->vertices);
      }
    }
  }

  /* ray/box intersection */
  __forceinline size_t intersectBox(const Ray8& ray, const avx3f& rdir, const avxf& ray_far, const BVH4MB::Node* node, const int i) 
  {
    if (unlikely(node->child[i]->isEmptyLeaf())) return 0;

    const avxf lower_x = avxf(node->lower_x[i]) + ray.time * avxf(node->lower_dx[i]);
    const avxf lower_y = avxf(node->lower_y[i]) + ray.time * avxf(node->lower_dy[i]);
    const avxf lower_z = avxf(node->lower_z[i]) + ray.time * avxf(node->lower_dz[i]);
    const avxf upper_x = avxf(node->upper_x[i]) + ray.time * avxf(node->upper_dx[i]);
    const avxf upper_y = avxf(node->upper_y[i]) + ray.time * avxf(node->upper_dy[i]);
    const avxf upper_z = avxf(node->upper_z[i]) + ray.time * avxf(node->upper_dz[i]);

    const avxf dminx = (lower_x - ray.org.x) * rdir.x;
    const avxf dminy = (lower_y - ray.org.y) * rdir.y;
    const avxf dminz = (lower_z - ray.org.z) * rdir.z;
    const avxf dmaxx = (upper_x - ray.org.x) * rdir.x;
    const avxf dmaxy = (upper_y - ray.org.y) * rdir.y;
    const avxf dmaxz = (upper_z - ray.org.z) * rdir.z;
    
    const avxf dlowerx = min(dminx,dmaxx);
    const avxf dlowery = min(dminy,dmaxy);
    const avxf dlowerz = min(dminz,dmaxz);
    
    const avxf dupperx = max(dminx,dmaxx);
    const avxf duppery = max(dminy,dmaxy);
    const avxf dupperz = max(dminz,dmaxz);
    
    const avxf near = max(dlowerx,dlowery,dlowerz,ray.tnear);
    const avxf far  = min(dupperx,duppery,dupperz,ray_far );
    
    return movemask(near <= far);
  }
  
  template<typename TriangleIntersector>
  __m256 BVH4MBIntersector8Chunk<TriangleIntersector>::occluded(const BVH4MBIntersector8Chunk* This, Ray8& ray, const __m256 valid_i)
  {
    avxb valid = valid_i;
    STAT3(shadow.travs,1,popcnt(valid),4);
    avxb terminated = !valid;
    const BVH4MB* bvh = This->bvh;
    
    BVH4MB::Base* stack[2+3*BVH4MB::maxDepth];
    BVH4MB::Base** stackPtr = stack+1; //!< current stack pointer
    stack[0] = bvh->root; 

    /* let terminated rays miss all boxes */
    avx3f rdir = rcp_safe(ray.dir);
    avxf rayFar  = select(terminated,avxf(neg_inf),ray.tfar);

    while (true)
    {
      /* finish when the stack is empty */
      if (unlikely(stackPtr == stack)) break;
      BVH4MB::Base* cur  = *(--stackPtr);
        
      /* this is an inner node */
      if (likely(cur->isNode()))
      {
        STAT3(shadow.trav_nodes,1,popcnt(valid),4);
        const BVH4MB::Node* node = cur->node();
        if (unlikely(intersectBox(ray,rdir,rayFar,node,0) != 0)) { *stackPtr = node->child[0]; stackPtr++; }
        if (unlikely(intersectBox(ray,rdir,rayFar,node,1) != 0)) { *stackPtr = node->child[1]; stackPtr++; }
        if (unlikely(intersectBox(ray,rdir,rayFar,node,2) != 0)) { *stackPtr = node->child[2]; stackPtr++; }
        if (unlikely(intersectBox(ray,rdir,rayFar,node,3) != 0)) { *stackPtr = node->child[3]; stackPtr++; }
      }

      /* this is a leaf node */
      else 
      {
        STAT3(shadow.trav_leaves,1,popcnt(valid),4);
        size_t num; Triangle* tri = (Triangle*) cur->leaf(num);
        for (size_t i=0; i<num; i++) {
          terminated |= TriangleIntersector::occluded(valid,ray,tri[i],bvh->vertices);
          if (all(terminated)) return terminated;
        }

        /* let terminated rays miss all boxes */
        rayFar = select(terminated,avxf(neg_inf),rayFar);
      }
    }
    return terminated;
  }

  void BVH4MBIntersector8ChunkRegister () 
  {
    //TriangleMesh::intersectors8.add("bvh4mb","triangle1i","chunk","moeller" ,false,BVH4MBIntersector8Chunk<Triangle1iIntersector8<Intersector8MoellerTrumbore> >::create);
    //TriangleMesh::intersectors8.add("bvh4mb","triangle1i","chunk","pluecker",true ,BVH4MBIntersector8Chunk<Triangle1iIntersector8<Intersector8Pluecker> >::create);
    //TriangleMesh::intersectors8.add("bvh4mb","triangle1v","chunk","moeller" ,false,BVH4MBIntersector8Chunk<Triangle1vIntersector8<Intersector8MoellerTrumbore> >::create);
    //TriangleMesh::intersectors8.add("bvh4mb","triangle1v","chunk","pluecker",true ,BVH4MBIntersector8Chunk<Triangle1vIntersector8<Intersector8Pluecker> >::create);
    TriangleMesh::intersectors8.add("bvh4mb","triangle4i","chunk","moeller" ,false,BVH4MBIntersector8Chunk<Triangle4iIntersector8<Intersector8MoellerTrumbore> >::create);
    TriangleMesh::intersectors8.add("bvh4mb","triangle4i","chunk","pluecker",true ,BVH4MBIntersector8Chunk<Triangle4iIntersector8<Intersector8Pluecker> >::create);
    //TriangleMesh::intersectors8.add("bvh4mb","triangle4v","chunk","moeller" ,false,BVH4MBIntersector8Chunk<Triangle4vIntersector8<Intersector8MoellerTrumbore> >::create);
    //TriangleMesh::intersectors8.add("bvh4mb","triangle4v","chunk","pluecker",true ,BVH4MBIntersector8Chunk<Triangle4vIntersector8<Intersector8Pluecker> >::create);
    //TriangleMesh::intersectors8.add("bvh4mb","triangle1" ,"chunk","moeller" ,true ,BVH4MBIntersector8Chunk<Triangle1Intersector8MoellerTrumbore>::create);
    //TriangleMesh::intersectors8.add("bvh4mb","triangle4" ,"chunk","moeller" ,true ,BVH4MBIntersector8Chunk<Triangle4Intersector8MoellerTrumbore>::create);
    //TriangleMesh::intersectors8.add("bvh4mb","triangle8" ,"chunk","moeller" ,true ,BVH4MBIntersector8Chunk<Triangle8Intersector8MoellerTrumbore>::create);
    //VirtualScene::intersectors8.add("bvh4mb","virtual"   ,"chunk","virtual" ,true ,BVH4MBIntersector8Chunk<VirtualObjectIntersector8>::create);
  }
}
#endif

