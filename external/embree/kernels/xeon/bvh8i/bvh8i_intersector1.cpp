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

#include "bvh8i_intersector1.h"
#include "../common/stack_item.h"
#include "bvh4i/bvh4i_builder_util.h"

#include "bvh8i.h"

#if defined(__AVX__)
#include "geometry/triangle8_intersector1_moeller.h"
#endif

namespace embree
{
  namespace isa
  {
// #if defined(__AVX2__)
// #if defined(__X86_64__)
//     __forceinline size_t bitscan(size_t mask) { return _tzcnt_u64(mask); }
// #else
//     __forceinline size_t bitscan(size_t mask) { return _tzcnt_u32(mask); }
// #endif
// #else
//     __forceinline size_t bitscan(size_t mask) { return __bsf(mask); }
// #endif



    
#define DBG(x) 
    
    template<typename TriangleIntersector>
    void BVH8iIntersector1<TriangleIntersector>::intersect(const BVH8i* bvh, Ray& ray)
    {
#if defined(__AVX__)
      /*! stack state */
      StackItemInt64 stack[1+3*BVH4i::maxDepth];  //!< stack of nodes 
      StackItemInt64* stackPtr = stack+1;        //!< current stack pointer
      stack[0].ptr  = bvh->root;
      stack[0].dist = neg_inf;
      
      /*! offsets to select the side that becomes the lower or upper bound */
      const size_t nearX = ray.dir.x >= 0.0f ? 0*sizeof(avxf) : 1*sizeof(avxf);
      const size_t nearY = ray.dir.y >= 0.0f ? 2*sizeof(avxf) : 3*sizeof(avxf);
      const size_t nearZ = ray.dir.z >= 0.0f ? 4*sizeof(avxf) : 5*sizeof(avxf);
      
      /*! load the ray into SIMD registers */
      const avx3f norg(-ray.org.x,-ray.org.y,-ray.org.z);
      const Vec3fa ray_rdir = rcp_safe(ray.dir);
      const avx3f rdir(ray_rdir.x,ray_rdir.y,ray_rdir.z);
      const Vec3fa ray_org_rdir = ray.org*ray_rdir;
      const avx3f org_rdir(ray_org_rdir.x,ray_org_rdir.y,ray_org_rdir.z);
      const avxf  rayNear(ray.tnear);
      avxf rayFar(ray.tfar);
      
      const void* nodePtr = bvh->nodePtr();
      const void* triPtr  = bvh->triPtr();

      /* pop loop */
      while (true) pop:
      {
        /*! pop next node */
        if (unlikely(stackPtr == stack)) break;
        stackPtr--;
        NodeRef cur = NodeRef(stackPtr->ptr);
        
        /*! if popped node is too far, pop next one */
        if (unlikely(*(float*)&stackPtr->dist > ray.tfar))
          continue;
        
        /* downtraversal loop */
        while (true)
        {
          /*! stop if we found a leaf */
          if (unlikely(cur.isLeaf())) break;
          STAT3(normal.trav_nodes,1,1,1);
          
          /*! single ray intersection with 4 boxes */
          const Node* node = (BVH8i::Node*)cur.node(nodePtr);
          const size_t farX  = nearX ^ sizeof(avxf), farY  = nearY ^ sizeof(avxf), farZ  = nearZ ^ sizeof(avxf);
          const avxf tNearX = msub(load8f((const char*)nodePtr+(size_t)cur+nearX), rdir.x, org_rdir.x);
          const avxf tNearY = msub(load8f((const char*)nodePtr+(size_t)cur+nearY), rdir.y, org_rdir.y);
          const avxf tNearZ = msub(load8f((const char*)nodePtr+(size_t)cur+nearZ), rdir.z, org_rdir.z);
          const avxf tFarX  = msub(load8f((const char*)nodePtr+(size_t)cur+farX ), rdir.x, org_rdir.x);
          const avxf tFarY  = msub(load8f((const char*)nodePtr+(size_t)cur+farY ), rdir.y, org_rdir.y);
          const avxf tFarZ  = msub(load8f((const char*)nodePtr+(size_t)cur+farZ ), rdir.z, org_rdir.z);
          
#if defined(__AVX2__)
          const avxf tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,rayNear));
          const avxf tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,rayFar ));
          const avxb vmask = cast(tNear) > cast(tFar);
          unsigned int mask = movemask(vmask)^0xff;
#else
          const avxf tNear = max(tNearX,tNearY,tNearZ,rayNear);
          const avxf tFar  = min(tFarX ,tFarY ,tFarZ ,rayFar);
          const avxb vmask = tNear <= tFar;
          size_t mask = movemask(vmask);
#endif
          
          /*! if no child is hit, pop next node */
          if (unlikely(mask == 0))
            goto pop;
          
          /*! one child is hit, continue with that child */
          size_t r = __bscf(mask);
          if (likely(mask == 0)) {
            cur = node->child(r);
            continue;
          }
          
          /*! two children are hit, push far child, and continue with closer child */
          NodeRef c0 = node->child(r); const unsigned int d0 = *(unsigned int*)&tNear[r];
          r = __bscf(mask);
          NodeRef c1 = node->child(r); const unsigned int d1 = *(unsigned int*)&tNear[r];
          if (likely(mask == 0)) {
            if (d0 < d1) { stackPtr->ptr = c1; stackPtr->dist = d1; stackPtr++; cur = c0; continue; }
            else         { stackPtr->ptr = c0; stackPtr->dist = d0; stackPtr++; cur = c1; continue; }
          }
          
          /*! Here starts the slow path for 3 or 4 hit children. We push
           *  all nodes onto the stack to sort them there. */
          stackPtr->ptr = c0; stackPtr->dist = d0; stackPtr++;
          stackPtr->ptr = c1; stackPtr->dist = d1; stackPtr++;
          
          /*! three children are hit, push all onto stack and sort 3 stack items, continue with closest child */
          r = __bscf(mask);
          NodeRef c = node->child(r); unsigned int d = *(unsigned int*)&tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
          if (likely(mask == 0)) {
            sort(stackPtr[-1],stackPtr[-2],stackPtr[-3]);
            cur = (NodeRef) stackPtr[-1].ptr; stackPtr--;
            continue;
          }
          
          /*! four children are hit, push all onto stack and sort 4 stack items, continue with closest child */
          r = __bscf(mask);
          c = node->child(r); d = *(unsigned int*)&tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
	  if (likely(mask == 0)) {
	    sort(stackPtr[-1],stackPtr[-2],stackPtr[-3],stackPtr[-4]);
	    cur = (NodeRef) stackPtr[-1].ptr; stackPtr--;
	    continue;
	  }

	  while(1)
	    {
	      r = __bscf(mask);
	      c = node->child(r); d = *(unsigned int*)&tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
	      if (unlikely(mask == 0)) break;
	    }

	    cur = (NodeRef) stackPtr[-1].ptr; stackPtr--;
        }
        
        /*! this is a leaf node */
        STAT3(normal.trav_leaves,1,1,1);
        size_t num; Triangle* tri = (Triangle*) cur.leaf(triPtr,num);
        TriangleIntersector::intersect(ray,tri,num,bvh->geometry);
        rayFar = ray.tfar;
      }
#endif
    }

    template<typename TriangleIntersector>    
    void BVH8iIntersector1<TriangleIntersector>::occluded(const BVH8i* bvh, Ray& ray)
    {
#if defined(__AVX__)
      /*! stack state */
      NodeRef stack[1+3*BVH4i::maxDepth];  //!< stack of nodes that still need to get traversed
      NodeRef* stackPtr = stack+1;        //!< current stack pointer
      stack[0]  = bvh->root;
      
      /*! offsets to select the side that becomes the lower or upper bound */
      const size_t nearX = ray.dir.x >= 0.0f ? 0*sizeof(avxf) : 1*sizeof(avxf);
      const size_t nearY = ray.dir.y >= 0.0f ? 2*sizeof(avxf) : 3*sizeof(avxf);
      const size_t nearZ = ray.dir.z >= 0.0f ? 4*sizeof(avxf) : 5*sizeof(avxf);
      
      /*! load the ray into SIMD registers */
      const avx3f norg(-ray.org.x,-ray.org.y,-ray.org.z);
      const Vec3fa ray_rdir = rcp_safe(ray.dir);
      const avx3f rdir(ray_rdir.x,ray_rdir.y,ray_rdir.z);
      const Vec3fa ray_org_rdir = ray.org*ray_rdir;
      const avx3f org_rdir(ray_org_rdir.x,ray_org_rdir.y,ray_org_rdir.z);
      const avxf  rayNear(ray.tnear);
      avxf rayFar(ray.tfar);

      const void* nodePtr = bvh->nodePtr();
      const void* triPtr  = bvh->triPtr();

      /* pop loop */
      while (true) pop:
      {
        /*! pop next node */
        if (unlikely(stackPtr == stack)) break;
        stackPtr--;
        NodeRef cur = (NodeRef) *stackPtr;
        
        /* downtraversal loop */
        while (true)
        {
          /*! stop if we found a leaf */
          if (unlikely(cur.isLeaf())) break;
          STAT3(shadow.trav_nodes,1,1,1);
          
          /*! single ray intersection with 4 boxes */
          const Node* node = (BVH8i::Node*)cur.node(nodePtr);
          const size_t farX  = nearX ^ sizeof(avxf), farY  = nearY ^ sizeof(avxf), farZ  = nearZ ^ sizeof(avxf);
          const avxf tNearX = msub(load8f((const char*)nodePtr+(size_t)cur+nearX), rdir.x, org_rdir.x);
          const avxf tNearY = msub(load8f((const char*)nodePtr+(size_t)cur+nearY), rdir.y, org_rdir.y);
          const avxf tNearZ = msub(load8f((const char*)nodePtr+(size_t)cur+nearZ), rdir.z, org_rdir.z);
          const avxf tFarX  = msub(load8f((const char*)nodePtr+(size_t)cur+farX ), rdir.x, org_rdir.x);
          const avxf tFarY  = msub(load8f((const char*)nodePtr+(size_t)cur+farY ), rdir.y, org_rdir.y);
          const avxf tFarZ  = msub(load8f((const char*)nodePtr+(size_t)cur+farZ ), rdir.z, org_rdir.z);
          
#if defined(__AVX2__)
          const avxf tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,rayNear));
          const avxf tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,rayFar ));
          const avxb vmask = cast(tNear) > cast(tFar);
          unsigned int mask = movemask(vmask)^0xff;
#else
          const avxf tNear = max(tNearX,tNearY,tNearZ,rayNear);
          const avxf tFar  = min(tFarX ,tFarY ,tFarZ ,rayFar);
          const avxb vmask = tNear <= tFar;
          size_t mask = movemask(vmask);
#endif
          
          /*! if no child is hit, pop next node */
          if (unlikely(mask == 0))
            goto pop;
          
          /*! one child is hit, continue with that child */
          size_t r = __bscf(mask);
          if (likely(mask == 0)) {
            cur = node->child(r);
            continue;
          }
          
          /*! two children are hit, push far child, and continue with closer child */
          NodeRef c0 = node->child(r); const unsigned int d0 = *(unsigned int*)&tNear[r];
          r = __bscf(mask);
          NodeRef c1 = node->child(r); const unsigned int d1 = *(unsigned int*)&tNear[r];

          if (likely(mask == 0)) {
            if (d0 < d1) { *stackPtr = c1; stackPtr++; cur = c0; continue; }
            else         { *stackPtr = c0; stackPtr++; cur = c1; continue; }
          }
          *stackPtr = c0; stackPtr++;
          *stackPtr = c1; stackPtr++;
          
          /*! three children are hit */
          r = __bscf(mask);
          cur = node->child(r); *stackPtr = cur; stackPtr++;
          if (likely(mask == 0)) {
            stackPtr--;
            continue;
          }

	  while(1)
	    {
	      r = __bscf(mask);
	      NodeRef c = node->child(r); *stackPtr = c; stackPtr++;
	      if (unlikely(mask == 0)) break;
	    }
	    cur = (NodeRef) stackPtr[-1]; stackPtr--;
          
        }
        
        /*! this is a leaf node */
        STAT3(shadow.trav_leaves,1,1,1);
        size_t num; Triangle* tri = (Triangle*) cur.leaf(triPtr,num);
        if (TriangleIntersector::occluded(ray,tri,num,bvh->geometry)) {
          ray.geomID = 0;
          break;
        }
      }      
#endif
      AVX_ZERO_UPPER();
    }

#if defined(__AVX__)    
    DEFINE_INTERSECTOR1(BVH8iTriangle8Intersector1Moeller,BVH8iIntersector1<Triangle8Intersector1MoellerTrumbore>);
#endif
  }
}
