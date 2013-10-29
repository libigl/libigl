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

#ifndef __EMBREE_BVH4_INTERSECTOR1_H__
#define __EMBREE_BVH4_INTERSECTOR1_H__

#include "bvh4.h"
#include "../include/intersector1.h"
#include "../common/stack_item.h"

namespace embree
{
  /*! BVH4 Traverser. Single ray traversal implementation for a Quad BVH. */
  template<typename TriangleIntersector>
    class BVH4Intersector1 : public Intersector1
  {
    /* shortcuts for frequently used types */
    typedef typename TriangleIntersector::Triangle Triangle;
    typedef typename BVH4::NodeRef NodeRef;
    typedef typename BVH4::Node Node;
    typedef StackItemT<size_t> StackItem;

  public:
    BVH4Intersector1 (const BVH4* bvh) 
      : Intersector1((intersectFunc)intersect,(occludedFunc)occluded), bvh(bvh) {}

    static Intersector1* create(const Accel* bvh) { 
      return new BVH4Intersector1((const BVH4*)bvh); 
    }

    static void intersect(const BVH4Intersector1* This, Ray& ray);
    static bool occluded (const BVH4Intersector1* This, Ray& ray);

    template<typename Ray, typename T>
      __forceinline static void intersect1(const BVH4* bvh, NodeRef root, size_t k, Ray& ray, const T& ray_rdir)
    {
      /*! stack state */
      NodeRef popCur = root;            //!< pre-popped top node from the stack
      float popDist = neg_inf;              //!< pre-popped distance of top node from the stack
      StackItem stack[1+3*BVH4::maxDepth];  //!< stack of nodes that still need to get traversed
      StackItem* stackPtr = stack+1;        //!< current stack pointer
      
      /*! offsets to select the side that becomes the lower or upper bound */
      const size_t nearX = ray.dir.x[k] >= 0 ? 0*sizeof(ssef_m) : 1*sizeof(ssef_m);
      const size_t nearY = ray.dir.y[k] >= 0 ? 2*sizeof(ssef_m) : 3*sizeof(ssef_m);
      const size_t nearZ = ray.dir.z[k] >= 0 ? 4*sizeof(ssef_m) : 5*sizeof(ssef_m);
      const size_t farX  = nearX ^ 16;
      const size_t farY  = nearY ^ 16;
      const size_t farZ  = nearZ ^ 16;
      
      /*! load the ray into SIMD registers */
      const sse3f org(ray.org.x[k],ray.org.y[k],ray.org.z[k]);
      const sse3f rdir(ray_rdir.x[k],ray_rdir.y[k],ray_rdir.z[k]);
      const ssef rayNear(ray.tnear[k]);
      ssef rayFar(ray.tfar[k]);
      
      while (true)
      {
        /*! pop next node */
        if (unlikely(stackPtr == stack)) break;
        stackPtr--;
        NodeRef cur = popCur;
        
        /*! if popped node is too far, pop next one */
        if (unlikely(popDist > ray.tfar[k])) {
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
          const Node* node = cur.node(bvh->nodePtr());
          const ssef tNearX = (ssef((const char*)node+nearX) - org.x) * rdir.x;
          const ssef tNearY = (ssef((const char*)node+nearY) - org.y) * rdir.y;
          const ssef tNearZ = (ssef((const char*)node+nearZ) - org.z) * rdir.z;
          const ssef tNear = max(tNearX,tNearY,tNearZ,rayNear);
          const ssef tFarX = (ssef((const char*)node+farX) - org.x) * rdir.x;
          const ssef tFarY = (ssef((const char*)node+farY) - org.y) * rdir.y;
          const ssef tFarZ = (ssef((const char*)node+farZ) - org.z) * rdir.z;
          popCur = (NodeRef) stackPtr[-1].ptr;  //!< pre-pop of topmost stack item
          popDist = stackPtr[-1].dist;              //!< pre-pop of distance of topmost stack item
          const ssef tFar = min(tFarX,tFarY,tFarZ,rayFar);
          size_t _hit = movemask(tNear <= tFar);
          
          /*! if no child is hit, pop next node */
          if (unlikely(_hit == 0))
            continue;
          
          /*! one child is hit, continue with that child */
          size_t r = __bsf(_hit); _hit = __btc(_hit,r);
          if (likely(_hit == 0)) {
            cur = node->child(r);
            goto next;
          }
          
          /*! two children are hit, push far child, and continue with closer child */
          NodeRef c0 = node->child(r); const float d0 = tNear[r];
          r = __bsf(_hit); _hit = __btc(_hit,r);
          NodeRef c1 = node->child(r); const float d1 = tNear[r];
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
          NodeRef c = node->child(r); float d = tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
          if (likely(_hit == 0)) {
            sort(stackPtr[-1],stackPtr[-2],stackPtr[-3]);
            cur = (NodeRef) stackPtr[-1].ptr; stackPtr--;
            goto next;
          }
          
          /*! four children are hit, push all onto stack and sort 4 stack items, continue with closest child */
          r = __bsf(_hit); _hit = __btc(_hit,r);
          c = node->child(r); d = tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
          sort(stackPtr[-1],stackPtr[-2],stackPtr[-3],stackPtr[-4]);
          cur = (NodeRef) stackPtr[-1].ptr; stackPtr--;
          goto next;
        }
        
        /*! this is a leaf node */
        else 
        {
          STAT3(normal.trav_leaves,1,1,1);
          size_t num; Triangle* tri = (Triangle*) cur.leaf(bvh->triPtr(),num);
          for (size_t i=0; i<num; i++)
            TriangleIntersector::intersect(k,ray,tri[i],bvh->vertices);
          
          popCur = (NodeRef) stackPtr[-1].ptr;    //!< pre-pop of topmost stack item
          popDist = stackPtr[-1].dist;                //!< pre-pop of distance of topmost stack item
          rayFar = ray.tfar[k];
        }
      }
    }
  
    template<typename Ray, typename T>
      __forceinline static bool occluded1(const BVH4* bvh, NodeRef root, size_t k, Ray& ray, const T& ray_rdir)
    {
      /*! stack state */
      NodeRef stack[1+3*BVH4::maxDepth];  //!< stack of nodes that still need to get traversed
      NodeRef* stackPtr = stack+1;        //!< current stack pointer
      stack[0] = root;                        //!< push first node onto stack
      
      /*! offsets to select the side that becomes the lower or upper bound */
      const size_t nearX = ray.dir.x[k] >= 0 ? 0*sizeof(ssef_m) : 1*sizeof(ssef_m);
      const size_t nearY = ray.dir.y[k] >= 0 ? 2*sizeof(ssef_m) : 3*sizeof(ssef_m);
      const size_t nearZ = ray.dir.z[k] >= 0 ? 4*sizeof(ssef_m) : 5*sizeof(ssef_m);
      const size_t farX  = nearX ^ 16;
      const size_t farY  = nearY ^ 16;
      const size_t farZ  = nearZ ^ 16;
      
      /*! load the ray into SIMD registers */
      //const sse3f norg(-ray.org.x[k],-ray.org.y[k],-ray.org.z[k]);
      const sse3f org(ray.org.x[k],ray.org.y[k],ray.org.z[k]);
      const sse3f rdir(ray_rdir.x[k],ray_rdir.y[k],ray_rdir.z[k]);
      const ssef rayNear(ray.tnear[k]);
      const ssef rayFar (ray.tfar[k]);
      
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
          const Node* node = cur.node(bvh->nodePtr());
          const ssef tNearX = (ssef((const char*)node+nearX) - org.x) * rdir.x;
          const ssef tNearY = (ssef((const char*)node+nearY) - org.y) * rdir.y;
          const ssef tNearZ = (ssef((const char*)node+nearZ) - org.z) * rdir.z;
          const ssef tNear = max(tNearX,tNearY,tNearZ,rayNear);
          const ssef tFarX = (ssef((const char*)node+farX) - org.x) * rdir.x;
          const ssef tFarY = (ssef((const char*)node+farY) - org.y) * rdir.y;
          const ssef tFarZ = (ssef((const char*)node+farZ) - org.z) * rdir.z;
          const ssef tFar = min(tFarX,tFarY,tFarZ,rayFar);
          size_t _hit = movemask(tNear <= tFar);

          /*! push hit nodes onto stack */
          if (likely(_hit == 0)) continue;
          size_t r = __bsf(_hit); _hit = __btc(_hit,r);
          *stackPtr = node->child(r); stackPtr++;
          if (likely(_hit == 0)) continue;
          r = __bsf(_hit); _hit = __btc(_hit,r);
          *stackPtr = node->child(r); stackPtr++;
          if (likely(_hit == 0)) continue;
          r = __bsf(_hit); _hit = __btc(_hit,r);
          *stackPtr = node->child(r); stackPtr++;
          if (likely(_hit == 0)) continue;
          r = __bsf(_hit); _hit = __btc(_hit,r);
          *stackPtr = node->child(r); stackPtr++;
        }
        
        /*! this is a leaf node */
        else 
        {
          STAT3(shadow.trav_leaves,1,1,1);
          size_t num; Triangle* tri = (Triangle*) cur.leaf(bvh->triPtr(),num);
          for (size_t i=0; i<num; i++)
            if (TriangleIntersector::occluded(k,ray,tri[i],bvh->vertices)) {
              AVX_ZERO_UPPER();
              return true;
            }
        }
      }
      AVX_ZERO_UPPER();
      return false;
    }

  private:
    const BVH4* bvh;
  };
}

#endif
