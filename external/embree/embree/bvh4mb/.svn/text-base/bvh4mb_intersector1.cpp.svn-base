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

#include "bvh4mb_intersector1.h"
#include "../common/stack_item.h"
#include "../geometry/triangles.h"

namespace embree
{
  template<typename TriangleIntersector>
  void BVH4MBIntersector1<TriangleIntersector>::intersect(const BVH4MBIntersector1* This, Ray& ray)
  {
    AVX_ZERO_UPPER();
    STAT3(normal.travs,1,1,1);

    /*! stack state */
    const BVH4MB* bvh = This->bvh;
    Base* popCur  = bvh->root;              //!< pre-popped top node from the stack
    float popDist = neg_inf;                //!< pre-popped distance of top node from the stack
    StackItem stack[1+3*BVH4MB::maxDepth];  //!< stack of nodes that still need to get traversed
    StackItem* stackPtr = stack+1;          //!< current stack pointer

    /*! offsets to select the side that becomes the lower or upper bound */
    const size_t nearX = ray.dir.x >= 0 ? 0*2*sizeof(ssef_m) : 1*2*sizeof(ssef_m);
    const size_t nearY = ray.dir.y >= 0 ? 2*2*sizeof(ssef_m) : 3*2*sizeof(ssef_m);
    const size_t nearZ = ray.dir.z >= 0 ? 4*2*sizeof(ssef_m) : 5*2*sizeof(ssef_m);
    const size_t farX  = nearX ^ 32;
    const size_t farY  = nearY ^ 32;
    const size_t farZ  = nearZ ^ 32;

    /*! load the ray into SIMD registers */
    const sse3f norg(-ray.org.x,-ray.org.y,-ray.org.z);
    const Vector3f ray_rdir = rcp_safe(ray.dir);
    const sse3f rdir(ray_rdir.x,ray_rdir.y,ray_rdir.z);
    const ssef rayNear(ray.tnear);
    ssef rayFar(ray.tfar);
    
    while (true)
    {
      /*! pop next node */
      if (unlikely(stackPtr == stack)) break;
      stackPtr--;
      Base* cur = popCur;
      
      /*! if popped node is too far, pop next one */
      if (unlikely(popDist > ray.tfar)) {
        popCur  = (Base*)stackPtr[-1].ptr;
        popDist = stackPtr[-1].dist;
        continue;
      }

    next:

      /*! we mostly go into the inner node case */
      if (likely(cur->isNode()))
      {
        STAT3(normal.trav_nodes,1,1,1);

        /*! single ray intersection with 4 boxes */
        const Node* node = cur->node();
        const ssef_m* pNearX = (const ssef_m*)((const char*)node+nearX);
        const ssef_m* pNearY = (const ssef_m*)((const char*)node+nearY);
        const ssef_m* pNearZ = (const ssef_m*)((const char*)node+nearZ);
        const ssef tNearX = (norg.x + ssef(pNearX[0]) + ray.time*pNearX[1]) * rdir.x;
        const ssef tNearY = (norg.y + ssef(pNearY[0]) + ray.time*pNearY[1]) * rdir.y;
        const ssef tNearZ = (norg.z + ssef(pNearZ[0]) + ray.time*pNearZ[1]) * rdir.z;
        const ssef tNear = max(tNearX,tNearY,tNearZ,rayNear);
        const ssef_m* pFarX = (const ssef_m*)((const char*)node+farX);
        const ssef_m* pFarY = (const ssef_m*)((const char*)node+farY);
        const ssef_m* pFarZ = (const ssef_m*)((const char*)node+farZ);
        const ssef tFarX = (norg.x + ssef(pFarX[0]) + ray.time*pFarX[1]) * rdir.x;
        const ssef tFarY = (norg.y + ssef(pFarY[0]) + ray.time*pFarY[1]) * rdir.y;
        const ssef tFarZ = (norg.z + ssef(pFarZ[0]) + ray.time*pFarZ[1]) * rdir.z;
        popCur = (Base*) stackPtr[-1].ptr;      //!< pre-pop of topmost stack item
        popDist = stackPtr[-1].dist;            //!< pre-pop of distance of topmost stack item
        const ssef tFar = min(tFarX,tFarY,tFarZ,rayFar);
	size_t _hit = movemask(tNear <= tFar);

        /*! if no child is hit, pop next node */
        if (unlikely(_hit == 0))
          continue;

        /*! one child is hit, continue with that child */
        size_t r = __bsf(_hit); _hit = __btc(_hit,r);
        if (likely(_hit == 0)) {
          cur = node->child[r];
          goto next;
        }

        /*! two children are hit, push far child, and continue with closer child */
        Base* c0 = node->child[r]; const float d0 = tNear[r];
        r = __bsf(_hit); _hit = __btc(_hit,r);
        Base* c1 = node->child[r]; const float d1 = tNear[r];
        if (likely(_hit == 0)) {
          if (d0 < d1) { stackPtr->ptr = c1; stackPtr->dist = d1; stackPtr++; cur = c0; goto next; }
          else         { stackPtr->ptr = c0; stackPtr->dist = d0; stackPtr++; cur = c1; goto next; }
        }

        /*! Here starts the slow path for 3 or 4 hit children. We push
         *  all nodes onto the stack to sort them there. */
        stackPtr->ptr = c0; stackPtr->dist = d0; stackPtr++;
        stackPtr->ptr = c1; stackPtr->dist = d1; stackPtr++;

        /*! three children are hit, push all onto stack and sort 3 stack items, continue with closest child */
        r = __bsf(_hit); _hit = __btc(_hit,r);
        Base* c = node->child[r]; float d = tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
        if (likely(_hit == 0)) {
          sort(stackPtr[-1],stackPtr[-2],stackPtr[-3]);
          cur = (Base*) stackPtr[-1].ptr; stackPtr--;
          goto next;
        }

        /*! four children are hit, push all onto stack and sort 4 stack items, continue with closest child */
        r = __bsf(_hit); _hit = __btc(_hit,r);
        c = node->child[r]; d = tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
        sort(stackPtr[-1],stackPtr[-2],stackPtr[-3],stackPtr[-4]);
        cur = (Base*) stackPtr[-1].ptr; stackPtr--;
        goto next;
      }

      /*! this is a leaf node */
      else 
      {
        STAT3(normal.trav_leaves,1,1,1);
        size_t num; Triangle* tri = (Triangle*) cur->leaf(num);
        for (size_t i=0; i<num; i++)
          TriangleIntersector::intersect(ray,tri[i],bvh->vertices);

        popCur = (Base*) stackPtr[-1].ptr;  //!< pre-pop of topmost stack item
        popDist = stackPtr[-1].dist;        //!< pre-pop of distance of topmost stack item
        rayFar = ray.tfar;
      }
    }
    AVX_ZERO_UPPER();
  }

  template<typename TriangleIntersector>
  bool BVH4MBIntersector1<TriangleIntersector>::occluded(const BVH4MBIntersector1* This, Ray& ray)
  {
    AVX_ZERO_UPPER();
    STAT3(shadow.travs,1,1,1);

    /*! stack state */
    const BVH4MB* bvh = This->bvh;
    Base* stack[1+3*BVH4MB::maxDepth];  //!< stack of nodes that still need to get traversed
    Base** stackPtr = stack+1;          //!< current stack pointer
    stack[0] = bvh->root;               //!< push first node onto stack

    /*! offsets to select the side that becomes the lower or upper bound */
    const size_t nearX = (ray.dir.x >= 0) ? 0*2*sizeof(ssef_m) : 1*2*sizeof(ssef_m);
    const size_t nearY = (ray.dir.y >= 0) ? 2*2*sizeof(ssef_m) : 3*2*sizeof(ssef_m);
    const size_t nearZ = (ray.dir.z >= 0) ? 4*2*sizeof(ssef_m) : 5*2*sizeof(ssef_m);
    const size_t farX  = nearX ^ 32;
    const size_t farY  = nearY ^ 32;
    const size_t farZ  = nearZ ^ 32;

    /*! load the ray into SIMD registers */
    const sse3f norg(-ray.org.x,-ray.org.y,-ray.org.z);
    const Vector3f ray_rdir = rcp_safe(ray.dir);
    const sse3f rdir(ray_rdir.x,ray_rdir.y,ray_rdir.z);
    const ssef rayNear(ray.tnear);
    const ssef rayFar (ray.tfar);
    
    /*! pop node from stack */
    while (true)
    {
      /* finish when the stack is empty */
      if (unlikely(stackPtr == stack)) break;
      Base* cur = *(--stackPtr);

      /*! this is an inner node */
      if (likely(cur->isNode()))
      {
        STAT3(shadow.trav_nodes,1,1,1);
        
        /*! single ray intersection with 4 boxes */
        const Node* node = cur->node();
        const ssef_m* pNearX = (const ssef_m*)((const char*)node+nearX);
        const ssef_m* pNearY = (const ssef_m*)((const char*)node+nearY);
        const ssef_m* pNearZ = (const ssef_m*)((const char*)node+nearZ);
        const ssef tNearX = (norg.x + ssef(pNearX[0]) + ray.time*pNearX[1]) * rdir.x;
        const ssef tNearY = (norg.y + ssef(pNearY[0]) + ray.time*pNearY[1]) * rdir.y;
        const ssef tNearZ = (norg.z + ssef(pNearZ[0]) + ray.time*pNearZ[1]) * rdir.z;
        const ssef tNear = max(tNearX,tNearY,tNearZ,rayNear);
        const ssef_m* pFarX = (const ssef_m*)((const char*)node+farX);
        const ssef_m* pFarY = (const ssef_m*)((const char*)node+farY);
        const ssef_m* pFarZ = (const ssef_m*)((const char*)node+farZ);
        const ssef tFarX = (norg.x + ssef(pFarX[0]) + ray.time*pFarX[1]) * rdir.x;
        const ssef tFarY = (norg.y + ssef(pFarY[0]) + ray.time*pFarY[1]) * rdir.y;
        const ssef tFarZ = (norg.z + ssef(pFarZ[0]) + ray.time*pFarZ[1]) * rdir.z;
        const ssef tFar = min(tFarX,tFarY,tFarZ,rayFar);
	size_t _hit = movemask(tNear <= tFar);

        /*! push hit nodes onto stack */
        if (likely(_hit == 0)) continue;
        size_t r = __bsf(_hit); _hit = __btc(_hit,r);
        *stackPtr = node->child[r]; stackPtr++;
        if (likely(_hit == 0)) continue;
        r = __bsf(_hit); _hit = __btc(_hit,r);
        *stackPtr = node->child[r]; stackPtr++;
        if (likely(_hit == 0)) continue;
        r = __bsf(_hit); _hit = __btc(_hit,r);
        *stackPtr = node->child[r]; stackPtr++;
        if (likely(_hit == 0)) continue;
        r = __bsf(_hit); _hit = __btc(_hit,r);
        *stackPtr = node->child[r]; stackPtr++;
      }

      /*! this is a leaf node */
      else 
      {
        STAT3(shadow.trav_leaves,1,1,1);
        size_t num; Triangle* tri = (Triangle*) cur->leaf(num);
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

  void BVH4MBIntersector1Register () 
  {
    //TriangleMesh::intersectors1.add("bvh4mb","triangle1i","fast","moeller" ,false,BVH4MBIntersector1<Triangle1iIntersector1MB<Intersector1MoellerTrumbore> >::create);
    //TriangleMesh::intersectors1.add("bvh4mb","triangle1i","fast","pluecker",true ,BVH4MBIntersector1<Triangle1iIntersector1MB<Intersector1Pluecker> >::create);
    TriangleMesh::intersectors1.add("bvh4mb","triangle4i","fast","moeller" ,false,BVH4MBIntersector1<Triangle4iIntersector1MB<Intersector1MoellerTrumbore> >::create);
    TriangleMesh::intersectors1.add("bvh4mb","triangle4i","fast","pluecker",true ,BVH4MBIntersector1<Triangle4iIntersector1MB<Intersector1Pluecker> >::create);
    TriangleMesh::intersectors1.setAccelDefaultTraverser("bvh4mb","fast");

    //VirtualScene::intersectors1.add("bvh4mb","virtual"   ,"fast","virtual" ,true ,BVH4MBIntersector1<VirtualObjectIntersector1>::create);
    //VirtualScene::intersectors1.setAccelDefaultTraverser("bvh4mb","fast");
  }
}
