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

#include "sys/platform.h"

#if defined(__AVX__)

#include "bvh8_intersector1.h"
#include "../common/stack_item.h"
#include "../geometry/triangles.h"

namespace embree
{
  template<typename TriangleIntersector>
  void BVH8Intersector1<TriangleIntersector>::intersect(const BVH8Intersector1* This, Ray& ray)
  {
    AVX_ZERO_UPPER();
    STAT3(normal.travs,1,1,1);
    
    /*! stack state */
    const BVH8* bvh = This->bvh;
    NodeRef popCur = bvh->root;       //!< pre-popped top node from the stack
    float popDist = neg_inf;              //!< pre-popped distance of top node from the stack
    StackItem stack[1+7*BVH8::maxDepth];  //!< stack of nodes that still need to get traversed
    StackItem* stackPtr = stack+1;        //!< current stack pointer

    /*! offsets to select the side that becomes the lower or upper bound */
    const size_t nearX = ray.dir.x >= 0 ? 0*sizeof(avxf_m) : 1*sizeof(avxf_m);
    const size_t nearY = ray.dir.y >= 0 ? 2*sizeof(avxf_m) : 3*sizeof(avxf_m);
    const size_t nearZ = ray.dir.z >= 0 ? 4*sizeof(avxf_m) : 5*sizeof(avxf_m);
    const size_t farX  = nearX ^ 32;
    const size_t farY  = nearY ^ 32;
    const size_t farZ  = nearZ ^ 32;

    /*! load the ray into SIMD registers */
    const Vector3f ray_rdir = rcp_safe(ray.dir);
    const avx3f norg(-ray.org.x,-ray.org.y,-ray.org.z);
    const avx3f rdir(ray_rdir.x,ray_rdir.y,ray_rdir.z);
    const avxf rayNear(ray.tnear);
    avxf rayFar(ray.tfar);
    
    while (true)
    {
      /*! pop next node */
      if (unlikely(stackPtr == stack)) break;
      stackPtr--;
      NodeRef cur = popCur;
      
      /*! if popped node is too far, pop next one */
      if (unlikely(popDist > ray.tfar)) {
        popCur  = (NodeRef)stackPtr[-1].ptr;
        popDist = stackPtr[-1].dist;
        continue;
      }

    next:

      /*! we mostly go into the inner node case */
      if (likely(cur.isNode()))
      {
        STAT3(normal.trav_nodes,1,1,1);
	
        /*! single ray intersection with 4 boxes */
        const Node* node = cur.node((void*)NULL);
        const avxf tNearX = (norg.x + avxf((const char*)node+nearX)) * rdir.x;
        const avxf tNearY = (norg.y + avxf((const char*)node+nearY)) * rdir.y;
        const avxf tNearZ = (norg.z + avxf((const char*)node+nearZ)) * rdir.z;
        const avxf tNear = max(tNearX,tNearY,tNearZ,rayNear);
        const avxf tFarX = (norg.x + avxf((const char*)node+farX)) * rdir.x;
        const avxf tFarY = (norg.y + avxf((const char*)node+farY)) * rdir.y;
        const avxf tFarZ = (norg.z + avxf((const char*)node+farZ)) * rdir.z;
        popCur = (NodeRef) stackPtr[-1].ptr;  //!< pre-pop of topmost stack item
        popDist = stackPtr[-1].dist;              //!< pre-pop of distance of topmost stack item
        const avxf tFar = min(tFarX,tFarY,tFarZ,rayFar);
	size_t mask = movemask(tNear <= tFar);

        /*! if no child is hit, pop next node */
        if (unlikely(mask == 0))
          continue;

        /*! one child is hit, continue with that child */
        size_t r = __bsf(mask); mask = __btc(mask,r);
        if (likely(mask == 0)) {
          cur = node->child(r);
          goto next;
        }

        /*! two children are hit, push far child, and continue with closer child */
        NodeRef c0 = node->child(r); const float d0 = tNear[r];
        r = __bsf(mask); mask = __btc(mask,r);
        NodeRef c1 = node->child(r); const float d1 = tNear[r];
        if (likely(mask == 0)) {
          if (d0 < d1) { stackPtr->ptr = c1; stackPtr->dist = d1; stackPtr++; cur = c0; goto next; }
          else         { stackPtr->ptr = c0; stackPtr->dist = d0; stackPtr++; cur = c1; goto next; }
        }

        /*! Here starts the slow path for 3 or 4 hit children. We push
         *  all nodes onto the stack to sort them there. */
        stackPtr->ptr = c0; stackPtr->dist = d0; stackPtr++;
        stackPtr->ptr = c1; stackPtr->dist = d1; stackPtr++;

        /*! three children are hit, push all onto stack and sort 3 stack items, continue with closest child */
        r = __bsf(mask); mask = __btc(mask,r);
        NodeRef c = node->child(r); float d = tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
        if (likely(mask == 0)) {
          sort(stackPtr[-1],stackPtr[-2],stackPtr[-3]);
          cur = (NodeRef) stackPtr[-1].ptr; stackPtr--;
          goto next;
        }

        /*! four children are hit, push all onto stack and sort 4 stack items, continue with closest child */
        r = __bsf(mask); mask = __btc(mask,r);
        c = node->child(r); d = tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
        if (likely(mask == 0)) {
          sort(stackPtr[-1],stackPtr[-2],stackPtr[-3],stackPtr[-4]);
          cur = (NodeRef) stackPtr[-1].ptr; stackPtr--;
          goto next;
        }

        /*! Here starts the slow path for 5 to 8 hit children. We push
         *  all nodes onto the stack and sort with a loop. */
        StackItem* basePtr = stackPtr-4;
        r = __bsf(mask); mask = __btc(mask,r);
        c = node->child(r); d = tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
        if (__builtin_expect(mask == 0, true)) goto sort_full;

        r = __bsf(mask); mask = __btc(mask,r);
        c = node->child(r); d = tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
        if (__builtin_expect(mask == 0, true)) goto sort_full;

        r = __bsf(mask); mask = __btc(mask,r);
        c = node->child(r); d = tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
        if (__builtin_expect(mask == 0, true)) goto sort_full;

        r = __bsf(mask); mask = __btc(mask,r);
        c = node->child(r); d = tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
       
      sort_full: // 1.8%
        size_t N = stackPtr-basePtr;
        for (ssize_t j=N-1; j>0; j--)
          for (ssize_t i=0; i<j; i++)
            sort(basePtr[i],basePtr[i+1]);

        cur = (NodeRef) stackPtr[-1].ptr; stackPtr--; 
        goto next;
      }

      /*! this is a leaf node */
      else 
      {
        STAT3(normal.trav_leaves,1,1,1);
        size_t num; Triangle* tri = (Triangle*) cur.leaf(NULL,num);
        for (size_t i=0; i<num; i++)
          TriangleIntersector::intersect(ray,tri[i],bvh->vertices);

        popCur = (NodeRef) stackPtr[-1].ptr;    //!< pre-pop of topmost stack item
        popDist = stackPtr[-1].dist;                //!< pre-pop of distance of topmost stack item
        rayFar = ray.tfar;
      }
    }
    AVX_ZERO_UPPER();
  }

  template<typename TriangleIntersector>
  bool BVH8Intersector1<TriangleIntersector>::occluded(const BVH8Intersector1* This, Ray& ray)
  {
    AVX_ZERO_UPPER();
    STAT3(shadow.travs,1,1,1);

    /*! stack state */
    const BVH8* bvh = This->bvh;
    NodeRef stack[1+3*BVH8::maxDepth];  //!< stack of nodes that still need to get traversed
    NodeRef* stackPtr = stack+1;        //!< current stack pointer
    stack[0] = bvh->root;                   //!< push first node onto stack

    /*! offsets to select the side that becomes the lower or upper bound */
    const size_t nearX = (ray.dir.x >= 0) ? 0*sizeof(avxf_m) : 1*sizeof(avxf_m);
    const size_t nearY = (ray.dir.y >= 0) ? 2*sizeof(avxf_m) : 3*sizeof(avxf_m);
    const size_t nearZ = (ray.dir.z >= 0) ? 4*sizeof(avxf_m) : 5*sizeof(avxf_m);
    const size_t farX  = nearX ^ 32;
    const size_t farY  = nearY ^ 32;
    const size_t farZ  = nearZ ^ 32;

    /*! load the ray into SIMD registers */
    const Vector3f ray_rdir = rcp_safe(ray.dir);
    const avx3f norg(-ray.org.x,-ray.org.y,-ray.org.z);
    const avx3f rdir(ray_rdir.x,ray_rdir.y,ray_rdir.z);
    const avxf rayNear(ray.tnear);
    const avxf rayFar (ray.tfar);
    
    /*! pop node from stack */
    while (true)
    {
      /* finish when the stack is empty */
      if (unlikely(stackPtr == stack)) break;
      NodeRef cur = *(--stackPtr);

      /*! this is an inner node */
      if (likely(cur.isNode()))
      {
        STAT3(shadow.trav_nodes,1,1,1);
              
        /*! single ray intersection with 4 boxes */
        const Node* node = cur.node((void*)NULL);
        const avxf tNearX = (norg.x + avxf((const char*)node+nearX)) * rdir.x;
        const avxf tNearY = (norg.y + avxf((const char*)node+nearY)) * rdir.y;
        const avxf tNearZ = (norg.z + avxf((const char*)node+nearZ)) * rdir.z;
        const avxf tNear = max(tNearX,tNearY,tNearZ,rayNear);
        const avxf tFarX = (norg.x + avxf((const char*)node+farX)) * rdir.x;
        const avxf tFarY = (norg.y + avxf((const char*)node+farY)) * rdir.y;
        const avxf tFarZ = (norg.z + avxf((const char*)node+farZ)) * rdir.z;
        const avxf tFar = min(tFarX,tFarY,tFarZ,rayFar);
	size_t mask = movemask(tNear <= tFar);

        /*! push hit nodes onto stack */
        if (likely(mask == 0)) continue;
        size_t r = __bsf(mask); mask = __btc(mask,r);
        *stackPtr = node->child(r); stackPtr++;
        if (likely(mask == 0)) continue;
        r = __bsf(mask); mask = __btc(mask,r);
        *stackPtr = node->child(r); stackPtr++;
        if (likely(mask == 0)) continue;
        r = __bsf(mask); mask = __btc(mask,r);
        *stackPtr = node->child(r); stackPtr++;
        if (likely(mask == 0)) continue;
        r = __bsf(mask); mask = __btc(mask,r);
        *stackPtr = node->child(r); stackPtr++;
        if (likely(mask == 0)) continue;
        r = __bsf(mask); mask = __btc(mask,r);
        *stackPtr = node->child(r); stackPtr++;
        if (likely(mask == 0)) continue;
        r = __bsf(mask); mask = __btc(mask,r);
        *stackPtr = node->child(r); stackPtr++;
        if (likely(mask == 0)) continue;
        r = __bsf(mask); mask = __btc(mask,r);
        *stackPtr = node->child(r); stackPtr++;
        if (likely(mask == 0)) continue;
        r = __bsf(mask); mask = __btc(mask,r);
        *stackPtr = node->child(r); stackPtr++;
      }

      /*! this is a leaf node */
      else 
      {
        STAT3(shadow.trav_leaves,1,1,1);
        size_t num; Triangle* tri = (Triangle*) cur.leaf(NULL,num);
        for (size_t i=0; i<num; i++)
          if (TriangleIntersector::occluded(ray,tri[i],bvh->vertices)) {
            AVX_ZERO_UPPER();
            return true;
          }
      }
    }
    AVX_ZERO_UPPER();
    return false;
  }

  void BVH8Intersector1Register () 
  {
    TriangleMesh::intersectors1.add("bvh8","triangle1i","fast","moeller" ,false,BVH8Intersector1<Triangle1iIntersector1<Intersector1MoellerTrumbore> >::create);
    TriangleMesh::intersectors1.add("bvh8","triangle1i","fast","pluecker",true ,BVH8Intersector1<Triangle1iIntersector1<Intersector1Pluecker> >::create);
    TriangleMesh::intersectors1.add("bvh8","triangle1v","fast","moeller" ,false,BVH8Intersector1<Triangle1vIntersector1<Intersector1MoellerTrumbore> >::create);
    TriangleMesh::intersectors1.add("bvh8","triangle1v","fast","pluecker",true ,BVH8Intersector1<Triangle1vIntersector1<Intersector1Pluecker> >::create);
    TriangleMesh::intersectors1.add("bvh8","triangle4i","fast","moeller" ,false,BVH8Intersector1<Triangle4iIntersector1<Intersector1MoellerTrumbore> >::create);
    TriangleMesh::intersectors1.add("bvh8","triangle4i","fast","pluecker",true ,BVH8Intersector1<Triangle4iIntersector1<Intersector1Pluecker> >::create);
    TriangleMesh::intersectors1.add("bvh8","triangle4v","fast","moeller" ,false,BVH8Intersector1<Triangle4vIntersector1<Intersector1MoellerTrumbore> >::create);
    TriangleMesh::intersectors1.add("bvh8","triangle4v","fast","pluecker",true ,BVH8Intersector1<Triangle4vIntersector1<Intersector1Pluecker> >::create);
    TriangleMesh::intersectors1.add("bvh8","triangle1" ,"fast","moeller" ,true ,BVH8Intersector1<Triangle1Intersector1MoellerTrumbore>::create);
    TriangleMesh::intersectors1.add("bvh8","triangle4" ,"fast","moeller" ,true ,BVH8Intersector1<Triangle4Intersector1MoellerTrumbore>::create);
    TriangleMesh::intersectors1.add("bvh8","triangle8" ,"fast","moeller" ,true ,BVH8Intersector1<Triangle8Intersector1MoellerTrumbore>::create);
    TriangleMesh::intersectors1.setAccelDefaultTraverser("bvh8","fast");

    VirtualScene::intersectors1.add("bvh8","virtual"   ,"fast","virtual" ,true ,BVH8Intersector1<VirtualObjectIntersector1>::create);
    VirtualScene::intersectors1.setAccelDefaultTraverser("bvh8","fast");
  }
}

#endif
