// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
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

#include "bvh8_intersector1.h"
#include "geometry/triangle4_intersector1_moeller.h"
#include "geometry/triangle8_intersector1_moeller.h"

namespace embree
{ 
  namespace isa
  {
    template<typename PrimitiveIntersector>
    void BVH8Intersector1<PrimitiveIntersector>::intersect(const BVH8* bvh, Ray& ray)
    {
      /*! perform per ray precalculations required by the primitive intersector */
      Precalculations pre(ray);

      /*! stack state */
      StackItemInt32<NodeRef> stack[stackSize];  //!< stack of nodes 
      StackItemInt32<NodeRef>* stackPtr = stack+1;        //!< current stack pointer
      StackItemInt32<NodeRef>* stackEnd = stack+stackSize;
      stack[0].ptr = bvh->root;
      stack[0].dist = neg_inf;
      
      /*! load the ray into SIMD registers */
      const avx3f norg(-ray.org.x,-ray.org.y,-ray.org.z);
      const Vec3fa ray_rdir = rcp_safe(ray.dir);
      const avx3f rdir(ray_rdir.x,ray_rdir.y,ray_rdir.z);
      const Vec3fa ray_org_rdir = ray.org*ray_rdir;
      const avx3f org_rdir(ray_org_rdir.x,ray_org_rdir.y,ray_org_rdir.z);
      const avxf  ray_near(ray.tnear);
      avxf ray_far(ray.tfar);

      /*! offsets to select the side that becomes the lower or upper bound */
      const size_t nearX = ray_rdir.x >= 0.0f ? 0*sizeof(avxf) : 1*sizeof(avxf);
      const size_t nearY = ray_rdir.y >= 0.0f ? 2*sizeof(avxf) : 3*sizeof(avxf);
      const size_t nearZ = ray_rdir.z >= 0.0f ? 4*sizeof(avxf) : 5*sizeof(avxf);
      
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
          const size_t farX  = nearX ^ sizeof(avxf), farY  = nearY ^ sizeof(avxf), farZ  = nearZ ^ sizeof(avxf);
#if defined (__AVX2__)
          const avxf tNearX = msub(load8f((const char*)node+nearX), rdir.x, org_rdir.x);
          const avxf tNearY = msub(load8f((const char*)node+nearY), rdir.y, org_rdir.y);
          const avxf tNearZ = msub(load8f((const char*)node+nearZ), rdir.z, org_rdir.z);
          const avxf tFarX  = msub(load8f((const char*)node+farX ), rdir.x, org_rdir.x);
          const avxf tFarY  = msub(load8f((const char*)node+farY ), rdir.y, org_rdir.y);
          const avxf tFarZ  = msub(load8f((const char*)node+farZ ), rdir.z, org_rdir.z);
#else
          const avxf tNearX = (norg.x + load8f((const char*)node+nearX)) * rdir.x;
          const avxf tNearY = (norg.y + load8f((const char*)node+nearY)) * rdir.y;
          const avxf tNearZ = (norg.z + load8f((const char*)node+nearZ)) * rdir.z;
          const avxf tFarX  = (norg.x + load8f((const char*)node+farX )) * rdir.x;
          const avxf tFarY  = (norg.y + load8f((const char*)node+farY )) * rdir.y;
          const avxf tFarZ  = (norg.z + load8f((const char*)node+farZ )) * rdir.z;
#endif

#if defined(__AVX2__)
          const avxf tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,ray_near));
          const avxf tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,ray_far ));
          const avxb vmask = cast(tNear) > cast(tFar);
          size_t mask = movemask(vmask)^0xff;
#else
          const avxf tNear = max(tNearX,tNearY,tNearZ,ray_near);
          const avxf tFar  = min(tFarX ,tFarY ,tFarZ ,ray_far);
          const avxb vmask = tNear <= tFar;
          size_t mask = movemask(vmask);
#endif
          
          /*! if no child is hit, pop next node */
          if (unlikely(mask == 0))
            goto pop;
          
          /*! one child is hit, continue with that child */
          size_t r = __bscf(mask);
          if (likely(mask == 0)) {
            cur = node->child(r); cur.prefetch();
            assert(cur != BVH8::emptyNode);
            continue;
          }
          
          /*! two children are hit, push far child, and continue with closer child */
          NodeRef c0 = node->child(r); c0.prefetch(); const unsigned int d0 = ((unsigned int*)&tNear)[r];
          r = __bscf(mask);
          NodeRef c1 = node->child(r); c1.prefetch(); const unsigned int d1 = ((unsigned int*)&tNear)[r];
          assert(c0 != BVH8::emptyNode);
          assert(c1 != BVH8::emptyNode);
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
          NodeRef c = node->child(r); c.prefetch(); unsigned int d = ((unsigned int*)&tNear)[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
          assert(c != BVH8::emptyNode);
          if (likely(mask == 0)) {
            sort(stackPtr[-1],stackPtr[-2],stackPtr[-3]);
            cur = (NodeRef) stackPtr[-1].ptr; stackPtr--;
            continue;
          }
          
	  /*! four children are hit, push all onto stack and sort 4 stack items, continue with closest child */
          r = __bscf(mask);
          c = node->child(r); c.prefetch(); d = *(unsigned int*)&tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
	  if (likely(mask == 0)) {
	    sort(stackPtr[-1],stackPtr[-2],stackPtr[-3],stackPtr[-4]);
	    cur = (NodeRef) stackPtr[-1].ptr; stackPtr--;
	    continue;
	  }

	  /*! fallback case if more than 4 children are hit */
	  while (1)
	  {
	    r = __bscf(mask);
	    assert(stackPtr < stackEnd);
	    c = node->child(r); c.prefetch(); d = *(unsigned int*)&tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
	    if (unlikely(mask == 0)) break;
	  }
	  
	  cur = (NodeRef) stackPtr[-1].ptr; stackPtr--;
	}
        
        /*! this is a leaf node */
	assert(cur != BVH8::emptyNode);
        STAT3(normal.trav_leaves,1,1,1);
        size_t num; Primitive* prim = (Primitive*) cur.leaf(num);
        size_t lazy_node = 0;
        PrimitiveIntersector::intersect(pre,ray,prim,num,bvh->scene,lazy_node);
        ray_far = ray.tfar;

        if (unlikely(lazy_node)) {
          stackPtr->ptr = lazy_node;
          stackPtr->dist = inf;
          stackPtr++;
        }
      }

      AVX_ZERO_UPPER();
    }
    
    template<typename PrimitiveIntersector>
    void BVH8Intersector1<PrimitiveIntersector>::occluded(const BVH8* bvh, Ray& ray)
    {
      /*! perform per ray precalculations required by the primitive intersector */
      Precalculations pre(ray);

      /*! stack state */
      NodeRef stack[stackSize];  //!< stack of nodes that still need to get traversed
      NodeRef* stackPtr = stack+1;        //!< current stack pointer
      NodeRef* stackEnd = stack+stackSize;
      stack[0] = bvh->root;
            
      /*! load the ray into SIMD registers */
      const avx3f norg(-ray.org.x,-ray.org.y,-ray.org.z);
      const Vec3fa ray_rdir = rcp_safe(ray.dir);
      const avx3f rdir(ray_rdir.x,ray_rdir.y,ray_rdir.z);
      const Vec3fa ray_org_rdir = ray.org*ray_rdir;
      const avx3f org_rdir(ray_org_rdir.x,ray_org_rdir.y,ray_org_rdir.z);
      const avxf  ray_near(ray.tnear);
      avxf ray_far(ray.tfar);
      
      /*! offsets to select the side that becomes the lower or upper bound */
      const size_t nearX = ray_rdir.x >= 0 ? 0*sizeof(avxf) : 1*sizeof(avxf);
      const size_t nearY = ray_rdir.y >= 0 ? 2*sizeof(avxf) : 3*sizeof(avxf);
      const size_t nearZ = ray_rdir.z >= 0 ? 4*sizeof(avxf) : 5*sizeof(avxf);

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
          const size_t farX  = nearX ^ sizeof(avxf), farY  = nearY ^ sizeof(avxf), farZ  = nearZ ^ sizeof(avxf);
#if defined (__AVX2__)
          const avxf tNearX = msub(load8f((const char*)node+nearX), rdir.x, org_rdir.x);
          const avxf tNearY = msub(load8f((const char*)node+nearY), rdir.y, org_rdir.y);
          const avxf tNearZ = msub(load8f((const char*)node+nearZ), rdir.z, org_rdir.z);
          const avxf tFarX  = msub(load8f((const char*)node+farX ), rdir.x, org_rdir.x);
          const avxf tFarY  = msub(load8f((const char*)node+farY ), rdir.y, org_rdir.y);
          const avxf tFarZ  = msub(load8f((const char*)node+farZ ), rdir.z, org_rdir.z);
#else
          const avxf tNearX = (norg.x + load8f((const char*)node+nearX)) * rdir.x;
          const avxf tNearY = (norg.y + load8f((const char*)node+nearY)) * rdir.y;
          const avxf tNearZ = (norg.z + load8f((const char*)node+nearZ)) * rdir.z;
          const avxf tFarX  = (norg.x + load8f((const char*)node+farX )) * rdir.x;
          const avxf tFarY  = (norg.y + load8f((const char*)node+farY )) * rdir.y;
          const avxf tFarZ  = (norg.z + load8f((const char*)node+farZ )) * rdir.z;
#endif
          
#if defined(__AVX2__)
          const avxf tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,ray_near));
          const avxf tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,ray_far ));
          const avxb vmask = cast(tNear) > cast(tFar);
          size_t mask = movemask(vmask)^0xff;
#else
          const avxf tNear = max(tNearX,tNearY,tNearZ,ray_near);
          const avxf tFar  = min(tFarX ,tFarY ,tFarZ ,ray_far);
          const avxb vmask = tNear <= tFar;
          size_t mask = movemask(vmask);
#endif
          
          /*! if no child is hit, pop next node */
          if (unlikely(mask == 0))
            goto pop;
          
          /*! one child is hit, continue with that child */
          size_t r = __bscf(mask);
          if (likely(mask == 0)) {
            cur = node->child(r); cur.prefetch(); 
            assert(cur != BVH8::emptyNode);
            continue;
          }
          
          /*! two children are hit, push far child, and continue with closer child */
          NodeRef c0 = node->child(r); c0.prefetch(); const unsigned int d0 = ((unsigned int*)&tNear)[r];
          r = __bscf(mask);
          NodeRef c1 = node->child(r); c1.prefetch(); const unsigned int d1 = ((unsigned int*)&tNear)[r];
          assert(c0 != BVH8::emptyNode);
          assert(c1 != BVH8::emptyNode);
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
          cur = node->child(r); cur.prefetch(); *stackPtr = cur; stackPtr++;
          if (likely(mask == 0)) {
            stackPtr--;
            continue;
          }

	  /*! process more than three children */
	  while(1)
	  {
	    r = __bscf(mask);
	    NodeRef c = node->child(r); c.prefetch(); *stackPtr = c; stackPtr++;
	    if (unlikely(mask == 0)) break;
	  }
	  cur = (NodeRef) stackPtr[-1]; stackPtr--;
        }
        
        /*! this is a leaf node */
	assert(cur != BVH8::emptyNode);
        STAT3(shadow.trav_leaves,1,1,1);
        size_t num; Primitive* prim = (Primitive*) cur.leaf(num);
        size_t lazy_node = 0;
        if (PrimitiveIntersector::occluded(pre,ray,prim,num,bvh->scene,lazy_node)) {
          ray.geomID = 0;
          break;
        }

        if (unlikely(lazy_node)) {
          *stackPtr = (NodeRef)lazy_node;
          stackPtr++;
        }
      }
      AVX_ZERO_UPPER();
    }

    DEFINE_INTERSECTOR1(BVH8Triangle4Intersector1Moeller,BVH8Intersector1<LeafIterator1<Triangle4Intersector1MoellerTrumbore<LeafMode> > >);
    DEFINE_INTERSECTOR1(BVH8Triangle8Intersector1Moeller,BVH8Intersector1<LeafIterator1<Triangle8Intersector1MoellerTrumbore<LeafMode> > >);
  }
}
