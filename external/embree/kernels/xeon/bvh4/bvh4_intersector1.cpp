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

#include "bvh4_intersector1.h"

#include "geometry/triangle1_intersector1_moeller.h"
#include "geometry/triangle4_intersector1_moeller.h"
#if defined(__AVX__)
#include "geometry/triangle8_intersector1_moeller.h"
#endif
#include "geometry/triangle1v_intersector1_pluecker.h"
#include "geometry/triangle4v_intersector1_pluecker.h"
#include "geometry/triangle4i_intersector1.h"
#include "geometry/virtual_accel_intersector1.h"

namespace embree
{ 
  namespace isa
  {
    template<typename PrimitiveIntersector>
    void BVH4Intersector1<PrimitiveIntersector>::intersect(const BVH4* bvh, Ray& ray)
    {
      /*! stack state */
      StackItemInt32<NodeRef> stack[stackSize];  //!< stack of nodes 
      StackItemInt32<NodeRef>* stackPtr = stack+1;        //!< current stack pointer
      StackItemInt32<NodeRef>* stackEnd = stack+stackSize;
      stack[0].ptr = bvh->root;
      stack[0].dist = neg_inf;
      
      /*! offsets to select the side that becomes the lower or upper bound */
      const size_t nearX = ray.dir.x >= 0.0f ? 0*sizeof(ssef) : 1*sizeof(ssef);
      const size_t nearY = ray.dir.y >= 0.0f ? 2*sizeof(ssef) : 3*sizeof(ssef);
      const size_t nearZ = ray.dir.z >= 0.0f ? 4*sizeof(ssef) : 5*sizeof(ssef);
      
#if 0 // FIXME: why is this slower
      /*! load the ray */
      Vec3fa ray_org = ray.org;
      Vec3fa ray_dir = ray.dir;
      ssef ray_near  = max(ray.tnear,FLT_MIN); // we do not support negative tnear values in this kernel due to integer optimizations
      ssef ray_far   = ray.tfar; 
#if defined(__FIX_RAYS__)
      const float float_range = 0.1f*FLT_MAX;
      ray_org = clamp(ray_org,Vec3fa(-float_range),Vec3fa(+float_range));
      ray_dir = clamp(ray_dir,Vec3fa(-float_range),Vec3fa(+float_range));
      ray_far = min(ray_far,float(inf)); 
#endif
      const Vec3fa ray_rdir = rcp_safe(ray_dir);
      const sse3f org(ray_org), dir(ray_dir);
      const sse3f norg(-ray_org), rdir(ray_rdir), org_rdir(ray_org*ray_rdir);
#else
      /*! load the ray into SIMD registers */
      const sse3f norg(-ray.org.x,-ray.org.y,-ray.org.z);
      const Vec3fa ray_rdir = rcp_safe(ray.dir);
      const sse3f rdir(ray_rdir.x,ray_rdir.y,ray_rdir.z);
      const Vec3fa ray_org_rdir = ray.org*ray_rdir;
      const sse3f org_rdir(ray_org_rdir.x,ray_org_rdir.y,ray_org_rdir.z);
      const ssef  ray_near(ray.tnear);
      ssef ray_far(ray.tfar);
#endif

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
          const Node* node = cur.node();
          const size_t farX  = nearX ^ 16, farY  = nearY ^ 16, farZ  = nearZ ^ 16;
#if defined (__AVX2__)
          const ssef tNearX = msub(load4f((const char*)node+nearX), rdir.x, org_rdir.x);
          const ssef tNearY = msub(load4f((const char*)node+nearY), rdir.y, org_rdir.y);
          const ssef tNearZ = msub(load4f((const char*)node+nearZ), rdir.z, org_rdir.z);
          const ssef tFarX  = msub(load4f((const char*)node+farX ), rdir.x, org_rdir.x);
          const ssef tFarY  = msub(load4f((const char*)node+farY ), rdir.y, org_rdir.y);
          const ssef tFarZ  = msub(load4f((const char*)node+farZ ), rdir.z, org_rdir.z);
#else
          const ssef tNearX = (norg.x + load4f((const char*)node+nearX)) * rdir.x;
          const ssef tNearY = (norg.y + load4f((const char*)node+nearY)) * rdir.y;
          const ssef tNearZ = (norg.z + load4f((const char*)node+nearZ)) * rdir.z;
          const ssef tFarX  = (norg.x + load4f((const char*)node+farX )) * rdir.x;
          const ssef tFarY  = (norg.y + load4f((const char*)node+farY )) * rdir.y;
          const ssef tFarZ  = (norg.z + load4f((const char*)node+farZ )) * rdir.z;
#endif

#if defined(__SSE4_1__)
          const ssef tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,ray_near));
          const ssef tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,ray_far ));
          const sseb vmask = cast(tNear) > cast(tFar);
          size_t mask = movemask(vmask)^0xf;
#else
          const ssef tNear = max(tNearX,tNearY,tNearZ,ray_near);
          const ssef tFar  = min(tFarX ,tFarY ,tFarZ ,ray_far);
          const sseb vmask = tNear <= tFar;
          size_t mask = movemask(vmask);
#endif
          
          /*! if no child is hit, pop next node */
          if (unlikely(mask == 0))
            goto pop;
          
          /*! one child is hit, continue with that child */
          size_t r = __bscf(mask);
          if (likely(mask == 0)) {
            cur = node->child(r);
            assert(cur != BVH4::emptyNode);
            continue;
          }
          
          /*! two children are hit, push far child, and continue with closer child */
          NodeRef c0 = node->child(r); const unsigned int d0 = ((unsigned int*)&tNear)[r];
          r = __bscf(mask);
          NodeRef c1 = node->child(r); const unsigned int d1 = ((unsigned int*)&tNear)[r];
          assert(c0 != BVH4::emptyNode);
          assert(c1 != BVH4::emptyNode);
          if (likely(mask == 0)) {
            assert(stackPtr < stackEnd); 
            if (d0 < d1) { stackPtr->ptr = c1; stackPtr->dist = d1; stackPtr++; cur = c0; continue; }
            else         { stackPtr->ptr = c0; stackPtr->dist = d0; stackPtr++; cur = c1; continue; }
          }
          
          /*! Here starts the slow path for 3 or 4 hit children. We push
           *  all nodes onto the stack to sort them there. */
          assert(stackPtr < stackEnd); 
          stackPtr->ptr = c0; stackPtr->dist = d0; stackPtr++;
          assert(stackPtr < stackEnd); 
          stackPtr->ptr = c1; stackPtr->dist = d1; stackPtr++;
          
          /*! three children are hit, push all onto stack and sort 3 stack items, continue with closest child */
          assert(stackPtr < stackEnd); 
          r = __bscf(mask);
          NodeRef c = node->child(r); unsigned int d = ((unsigned int*)&tNear)[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
          assert(c != BVH4::emptyNode);
          if (likely(mask == 0)) {
            sort(stackPtr[-1],stackPtr[-2],stackPtr[-3]);
            cur = (NodeRef) stackPtr[-1].ptr; stackPtr--;
            continue;
          }
          
          /*! four children are hit, push all onto stack and sort 4 stack items, continue with closest child */
          assert(stackPtr < stackEnd); 
          r = __bscf(mask);
          c = node->child(r); d = *(unsigned int*)&tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
          assert(c != BVH4::emptyNode);
          sort(stackPtr[-1],stackPtr[-2],stackPtr[-3],stackPtr[-4]);
          cur = (NodeRef) stackPtr[-1].ptr; stackPtr--;
        }
        
        /*! this is a leaf node */
        STAT3(normal.trav_leaves,1,1,1);
        size_t num; Primitive* prim = (Primitive*) cur.leaf(num);
        PrimitiveIntersector::intersect(ray,prim,num,bvh->geometry);
        ray_far = ray.tfar;
      }
    }
    
    template<typename PrimitiveIntersector>
    void BVH4Intersector1<PrimitiveIntersector>::occluded(const BVH4* bvh, Ray& ray)
    {
      /*! stack state */
      NodeRef stack[stackSize];  //!< stack of nodes that still need to get traversed
      NodeRef* stackPtr = stack+1;        //!< current stack pointer
      NodeRef* stackEnd = stack+stackSize;
      stack[0] = bvh->root;
      
      /*! offsets to select the side that becomes the lower or upper bound */
      const size_t nearX = ray.dir.x >= 0 ? 0*sizeof(ssef) : 1*sizeof(ssef);
      const size_t nearY = ray.dir.y >= 0 ? 2*sizeof(ssef) : 3*sizeof(ssef);
      const size_t nearZ = ray.dir.z >= 0 ? 4*sizeof(ssef) : 5*sizeof(ssef);
      
#if 0 // FIXME: why is this slower
      /*! load the ray */
      Vec3fa ray_org = ray.org;
      Vec3fa ray_dir = ray.dir;
      ssef ray_near  = max(ray.tnear,FLT_MIN); // we do not support negative tnear values in this kernel due to integer optimizations
      ssef ray_far   = ray.tfar; 
#if defined(__FIX_RAYS__)
      const float float_range = 0.1f*FLT_MAX;
      ray_org = clamp(ray_org,Vec3fa(-float_range),Vec3fa(+float_range));
      ray_dir = clamp(ray_dir,Vec3fa(-float_range),Vec3fa(+float_range));
      ray_far = min(ray_far,float(inf)); 
#endif
      const Vec3fa ray_rdir = rcp_safe(ray_dir);
      const sse3f org(ray_org), dir(ray_dir);
      const sse3f norg(-ray_org), rdir(ray_rdir), org_rdir(ray_org*ray_rdir);
#else
      /*! load the ray into SIMD registers */
      const sse3f norg(-ray.org.x,-ray.org.y,-ray.org.z);
      const Vec3fa ray_rdir = rcp_safe(ray.dir);
      const sse3f rdir(ray_rdir.x,ray_rdir.y,ray_rdir.z);
      const Vec3fa ray_org_rdir = ray.org*ray_rdir;
      const sse3f org_rdir(ray_org_rdir.x,ray_org_rdir.y,ray_org_rdir.z);
      const ssef  ray_near(ray.tnear);
      ssef ray_far(ray.tfar);
#endif
      
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
          const Node* node = cur.node();
          const size_t farX  = nearX ^ 16, farY  = nearY ^ 16, farZ  = nearZ ^ 16;
#if defined (__AVX2__)
          const ssef tNearX = msub(load4f((const char*)node+nearX), rdir.x, org_rdir.x);
          const ssef tNearY = msub(load4f((const char*)node+nearY), rdir.y, org_rdir.y);
          const ssef tNearZ = msub(load4f((const char*)node+nearZ), rdir.z, org_rdir.z);
          const ssef tFarX  = msub(load4f((const char*)node+farX ), rdir.x, org_rdir.x);
          const ssef tFarY  = msub(load4f((const char*)node+farY ), rdir.y, org_rdir.y);
          const ssef tFarZ  = msub(load4f((const char*)node+farZ ), rdir.z, org_rdir.z);
#else
          const ssef tNearX = (norg.x + load4f((const char*)node+nearX)) * rdir.x;
          const ssef tNearY = (norg.y + load4f((const char*)node+nearY)) * rdir.y;
          const ssef tNearZ = (norg.z + load4f((const char*)node+nearZ)) * rdir.z;
          const ssef tFarX  = (norg.x + load4f((const char*)node+farX )) * rdir.x;
          const ssef tFarY  = (norg.y + load4f((const char*)node+farY )) * rdir.y;
          const ssef tFarZ  = (norg.z + load4f((const char*)node+farZ )) * rdir.z;
#endif
          
#if defined(__SSE4_1__)
          const ssef tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,ray_near));
          const ssef tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,ray_far ));
          const sseb vmask = cast(tNear) > cast(tFar);
          size_t mask = movemask(vmask)^0xf;
#else
          const ssef tNear = max(tNearX,tNearY,tNearZ,ray_near);
          const ssef tFar  = min(tFarX ,tFarY ,tFarZ ,ray_far);
          const sseb vmask = tNear <= tFar;
          size_t mask = movemask(vmask);
#endif
          
          /*! if no child is hit, pop next node */
          if (unlikely(mask == 0))
            goto pop;
          
          /*! one child is hit, continue with that child */
          size_t r = __bscf(mask);
          if (likely(mask == 0)) {
            cur = node->child(r);
            assert(cur != BVH4::emptyNode);
            continue;
          }
          
          /*! two children are hit, push far child, and continue with closer child */
          NodeRef c0 = node->child(r); const unsigned int d0 = ((unsigned int*)&tNear)[r];
          r = __bscf(mask);
          NodeRef c1 = node->child(r); const unsigned int d1 = ((unsigned int*)&tNear)[r];
          assert(c0 != BVH4::emptyNode);
          assert(c1 != BVH4::emptyNode);
          if (likely(mask == 0)) {
            assert(stackPtr < stackEnd);
            if (d0 < d1) { *stackPtr = c1; stackPtr++; cur = c0; continue; }
            else         { *stackPtr = c0; stackPtr++; cur = c1; continue; }
          }
          assert(stackPtr < stackEnd);
          *stackPtr = c0; stackPtr++;
          assert(stackPtr < stackEnd);
          *stackPtr = c1; stackPtr++;
          
          /*! three children are hit */
          r = __bscf(mask);
          cur = node->child(r); 
          assert(cur != BVH4::emptyNode);
          if (likely(mask == 0)) continue;
          assert(stackPtr < stackEnd);
          *stackPtr = cur; stackPtr++;
          
          /*! four children are hit */
          cur = node->child(3);
          assert(cur != BVH4::emptyNode);
        }
        
        /*! this is a leaf node */
        STAT3(shadow.trav_leaves,1,1,1);
        size_t num; Primitive* prim = (Primitive*) cur.leaf(num);
        if (PrimitiveIntersector::occluded(ray,prim,num,bvh->geometry)) {
          ray.geomID = 0;
          break;
        }
      }
      AVX_ZERO_UPPER();
    }

    DEFINE_INTERSECTOR1(BVH4Triangle1Intersector1Moeller,BVH4Intersector1<Triangle1Intersector1MoellerTrumbore>);
    DEFINE_INTERSECTOR1(BVH4Triangle4Intersector1Moeller,BVH4Intersector1<Triangle4Intersector1MoellerTrumbore>);
#if defined(__AVX__)
    DEFINE_INTERSECTOR1(BVH4Triangle8Intersector1Moeller,BVH4Intersector1<Triangle8Intersector1MoellerTrumbore>);
#endif
    DEFINE_INTERSECTOR1(BVH4Triangle1vIntersector1Pluecker,BVH4Intersector1<Triangle1vIntersector1Pluecker>);
    DEFINE_INTERSECTOR1(BVH4Triangle4vIntersector1Pluecker,BVH4Intersector1<Triangle4vIntersector1Pluecker>);
    DEFINE_INTERSECTOR1(BVH4Triangle4iIntersector1Pluecker,BVH4Intersector1<Triangle4iIntersector1Pluecker>);
    DEFINE_INTERSECTOR1(BVH4VirtualIntersector1,BVH4Intersector1<VirtualAccelIntersector1>);
  }
}
