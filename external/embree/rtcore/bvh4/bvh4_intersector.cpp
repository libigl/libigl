// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#include "bvh4_intersector.h"
#include "../common/stack_item.h"
#include "../triangle/triangles.h"

namespace embree
{
  template<typename TriangleIntersector>
  void BVH4Intersector<TriangleIntersector>::intersect(const Ray& ray, Hit& hit) const
  {
    AVX_ZERO_UPPER();
    STAT3(normal.travs,1,1,1);
    
    /*! stack state */
    BVH4::Base* popCur = bvh->root;       //!< pre-popped top node from the stack
    float popDist = neg_inf;              //!< pre-popped distance of top node from the stack
    StackItem stack[1+3*BVH4::maxDepth];  //!< stack of nodes that still need to get traversed
    StackItem* stackPtr = stack+1;        //!< current stack pointer

    /*! offsets to select the side that becomes the lower or upper bound */
    const size_t nearX = ray.dir.x >= 0 ? 0*sizeof(ssef) : 1*sizeof(ssef);
    const size_t nearY = ray.dir.y >= 0 ? 2*sizeof(ssef) : 3*sizeof(ssef);
    const size_t nearZ = ray.dir.z >= 0 ? 4*sizeof(ssef) : 5*sizeof(ssef);
    const size_t farX  = nearX ^ 16;
    const size_t farY  = nearY ^ 16;
    const size_t farZ  = nearZ ^ 16;

    /*! load the ray into SIMD registers */
    const sse3f norg(-ray.org.x,-ray.org.y,-ray.org.z);
    const sse3f rdir(ray.rdir.x,ray.rdir.y,ray.rdir.z);
    const ssef rayNear(ray.near);
    ssef rayFar(ray.far);
    hit.t = min(hit.t,ray.far);
    
    while (true)
    {
      /*! pop next node */
      if (unlikely(stackPtr == stack)) break;
      stackPtr--;
      BVH4::Base* cur = popCur;
      
      /*! if popped node is too far, pop next one */
      if (unlikely(popDist > hit.t)) {
        popCur  = (BVH4::Base*)stackPtr[-1].ptr;
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
        const ssef tNearX = (norg.x + *(ssef*)((const char*)node+nearX)) * rdir.x;
        const ssef tNearY = (norg.y + *(ssef*)((const char*)node+nearY)) * rdir.y;
        const ssef tNearZ = (norg.z + *(ssef*)((const char*)node+nearZ)) * rdir.z;
        const ssef tNear = max(tNearX,tNearY,tNearZ,rayNear);
        const ssef tFarX = (norg.x + *(ssef*)((const char*)node+farX)) * rdir.x;
        const ssef tFarY = (norg.y + *(ssef*)((const char*)node+farY)) * rdir.y;
        const ssef tFarZ = (norg.z + *(ssef*)((const char*)node+farZ)) * rdir.z;
        popCur = (BVH4::Base*) stackPtr[-1].ptr;  //!< pre-pop of topmost stack item
        popDist = stackPtr[-1].dist;              //!< pre-pop of distance of topmost stack item
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
        BVH4::Base* c0 = node->child[r]; const float d0 = tNear[r];
        r = __bsf(_hit); _hit = __btc(_hit,r);
        BVH4::Base* c1 = node->child[r]; const float d1 = tNear[r];
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
        BVH4::Base* c = node->child[r]; float d = tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
        if (likely(_hit == 0)) {
          sort(stackPtr[-1],stackPtr[-2],stackPtr[-3]);
          cur = (BVH4::Base*) stackPtr[-1].ptr; stackPtr--;
          goto next;
        }

        /*! four children are hit, push all onto stack and sort 4 stack items, continue with closest child */
        r = __bsf(_hit); _hit = __btc(_hit,r);
        c = node->child[r]; d = tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
        sort(stackPtr[-1],stackPtr[-2],stackPtr[-3],stackPtr[-4]);
        cur = (BVH4::Base*) stackPtr[-1].ptr; stackPtr--;
        goto next;
      }

      /*! this is a leaf node */
      else 
      {
        STAT3(normal.trav_leaves,1,1,1);
        size_t num; Triangle* tri = (Triangle*) cur->leaf(num);
        for (size_t i=0; i<num; i++)
          TriangleIntersector::intersect(ray,hit,tri[i],bvh->vertices);

        popCur = (BVH4::Base*) stackPtr[-1].ptr;    //!< pre-pop of topmost stack item
        popDist = stackPtr[-1].dist;                //!< pre-pop of distance of topmost stack item
        rayFar = hit.t;
      }
    }
    AVX_ZERO_UPPER();
  }

  template<typename TriangleIntersector>
  bool BVH4Intersector<TriangleIntersector>::occluded(const Ray& ray) const
  {
    AVX_ZERO_UPPER();
    STAT3(shadow.travs,1,1,1);

    /*! stack state */
    BVH4::Base* stack[1+3*BVH4::maxDepth];  //!< stack of nodes that still need to get traversed
    BVH4::Base** stackPtr = stack+1;        //!< current stack pointer
    stack[0] = bvh->root;                   //!< push first node onto stack

    /*! offsets to select the side that becomes the lower or upper bound */
    const size_t nearX = (ray.dir.x >= 0) ? 0*sizeof(ssef) : 1*sizeof(ssef);
    const size_t nearY = (ray.dir.y >= 0) ? 2*sizeof(ssef) : 3*sizeof(ssef);
    const size_t nearZ = (ray.dir.z >= 0) ? 4*sizeof(ssef) : 5*sizeof(ssef);
    const size_t farX  = nearX ^ 16;
    const size_t farY  = nearY ^ 16;
    const size_t farZ  = nearZ ^ 16;

    /*! load the ray into SIMD registers */
    const sse3f norg(-ray.org.x,-ray.org.y,-ray.org.z);
    const sse3f rdir(ray.rdir.x,ray.rdir.y,ray.rdir.z);
    const ssef rayNear(ray.near);
    const ssef rayFar (ray.far);
    
    /*! pop node from stack */
    while (true)
    {
      /* finish when the stack is empty */
      if (unlikely(stackPtr == stack)) break;
      BVH4::Base* cur = *(--stackPtr);

      /*! this is an inner node */
      if (likely(cur->isNode()))
      {
        STAT3(shadow.trav_nodes,1,1,1);
        
        /*! single ray intersection with 4 boxes */
        const Node* node = cur->node();
        const ssef tNearX = (norg.x + *(ssef*)((const char*)node+nearX)) * rdir.x;
        const ssef tNearY = (norg.y + *(ssef*)((const char*)node+nearY)) * rdir.y;
        const ssef tNearZ = (norg.z + *(ssef*)((const char*)node+nearZ)) * rdir.z;
        const ssef tNear = max(tNearX,tNearY,tNearZ,rayNear);
        const ssef tFarX = (norg.x + *(ssef*)((const char*)node+farX)) * rdir.x;
        const ssef tFarY = (norg.y + *(ssef*)((const char*)node+farY)) * rdir.y;
        const ssef tFarZ = (norg.z + *(ssef*)((const char*)node+farZ)) * rdir.z;
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

  /* explicit template instantiation */
  INSTANTIATE_TEMPLATE_BY_INTERSECTOR(BVH4Intersector);
}
