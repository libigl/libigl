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

#ifndef __EMBREE_BVH4_INTERSECTOR8_HYBRID_H__
#define __EMBREE_BVH4_INTERSECTOR8_HYBRID_H__

#include "bvh4.h"
#include "common/ray8.h"
#include "common/stack_item.h"

#define ENABLE_PREFETCHING

namespace embree
{
  namespace isa 
  {
    /*! BVH4 Traverser. Hybrid Packet traversal implementation for a Quad BVH. */
    template<typename PrimitiveIntersector8>
      class BVH4Intersector8Hybrid 
      {
	/* shortcuts for frequently used types */
	typedef typename PrimitiveIntersector8::Primitive Primitive;
	typedef typename BVH4::NodeRef NodeRef;
	typedef typename BVH4::Node Node;
	typedef StackItemT<NodeRef> StackItem;
	static const size_t stackSizeSingle = 1+3*BVH4::maxDepth;
	static const size_t stackSizeChunk = 4*BVH4::maxDepth+1;

      public:
	  static __forceinline void intersect1(const BVH4* bvh, NodeRef root, const size_t k, Ray8& ray, const avx3f &ray_org, const avx3f &ray_dir, const avx3f &ray_rdir, const avxf &ray_tnear, const avxf &ray_tfar, const avx3i& nearXYZ)
	  {
	    /*! stack state */
	    StackItemInt32<NodeRef> stack[stackSizeSingle];  //!< stack of nodes 
	    StackItemInt32<NodeRef>* stackPtr = stack + 1;        //!< current stack pointer
	    StackItemInt32<NodeRef>* stackEnd = stack + stackSizeSingle;
	    stack[0].ptr = root;
	    stack[0].dist = neg_inf;

	    /*! offsets to select the side that becomes the lower or upper bound */
	    const size_t nearX = nearXYZ.x[k];
	    const size_t nearY = nearXYZ.y[k];
	    const size_t nearZ = nearXYZ.z[k];

	    /*! load the ray into SIMD registers */
	    const sse3f org(ray_org.x[k], ray_org.y[k], ray_org.z[k]);
	    const sse3f rdir(ray_rdir.x[k], ray_rdir.y[k], ray_rdir.z[k]);
	    const sse3f org_rdir(org*rdir);
	    ssef rayNear(ray_tnear[k]), rayFar(ray_tfar[k]);

	    /* pop loop */
	    while (true) pop:
	      {
		/*! pop next node */
		if (unlikely(stackPtr == stack)) break;
		stackPtr--;
		NodeRef cur = NodeRef(stackPtr->ptr);

		/*! if popped node is too far, pop next one */
		if (unlikely(*(float*)&stackPtr->dist > ray.tfar[k]))
		  continue;

		/* downtraversal loop */
		while (true)
		  {
		    /*! stop if we found a leaf */
		    if (unlikely(cur.isLeaf())) break;
		    STAT3(normal.trav_nodes, 1, 1, 1);

		    /*! single ray intersection with 4 boxes */
		    const Node* node = cur.node();
		    const size_t farX = nearX ^ 16, farY = nearY ^ 16, farZ = nearZ ^ 16;
#if defined (__AVX2__)
		    const ssef tNearX = msub(load4f((const char*)node + nearX), rdir.x, org_rdir.x);
		    const ssef tNearY = msub(load4f((const char*)node + nearY), rdir.y, org_rdir.y);
		    const ssef tNearZ = msub(load4f((const char*)node + nearZ), rdir.z, org_rdir.z);
		    const ssef tFarX = msub(load4f((const char*)node + farX), rdir.x, org_rdir.x);
		    const ssef tFarY = msub(load4f((const char*)node + farY), rdir.y, org_rdir.y);
		    const ssef tFarZ = msub(load4f((const char*)node + farZ), rdir.z, org_rdir.z);
#else
		    const ssef tNearX = (load4f((const char*)node + nearX) - org.x) * rdir.x;
		    const ssef tNearY = (load4f((const char*)node + nearY) - org.y) * rdir.y;
		    const ssef tNearZ = (load4f((const char*)node + nearZ) - org.z) * rdir.z;
		    const ssef tFarX = (load4f((const char*)node + farX) - org.x) * rdir.x;
		    const ssef tFarY = (load4f((const char*)node + farY) - org.y) * rdir.y;
		    const ssef tFarZ = (load4f((const char*)node + farZ) - org.z) * rdir.z;
#endif

#if defined(__SSE4_1__)
		    const ssef tNear = maxi(maxi(tNearX, tNearY), maxi(tNearZ, rayNear));
		    const ssef tFar = mini(mini(tFarX, tFarY), mini(tFarZ, rayFar));
		    const sseb vmask = cast(tNear) > cast(tFar);
		    size_t mask = movemask(vmask) ^ 0xf;
#else
		    const ssef tNear = max(tNearX, tNearY, tNearZ, rayNear);
		    const ssef tFar = min(tFarX, tFarY, tFarZ, rayFar);
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

#if defined(__AVX2__) && defined(ENABLE_PREFETCHING)
		    prefetchL1(((char*)c0.node()) + 0);
		    prefetchL1(((char*)c0.node()) + 64);
		    prefetchL1(((char*)c1.node()) + 0);
		    prefetchL1(((char*)c1.node()) + 64);
#endif

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


		    assert(c0 != BVH4::emptyNode);
		    if (likely(mask == 0)) {
		      sort(stackPtr[-1], stackPtr[-2], stackPtr[-3]);
		      cur = (NodeRef)stackPtr[-1].ptr; stackPtr--;
		      continue;
		    }

		    /*! four children are hit, push all onto stack and sort 4 stack items, continue with closest child */
		    assert(stackPtr < stackEnd);
		    r = __bscf(mask);
		    c = node->child(r); d = ((unsigned int*)&tNear)[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
		    assert(c != BVH4::emptyNode);
		    sort(stackPtr[-1], stackPtr[-2], stackPtr[-3], stackPtr[-4]);
		    cur = (NodeRef)stackPtr[-1].ptr; stackPtr--;
		  }

		/*! this is a leaf node */
		STAT3(normal.trav_leaves, 1, 1, 1);
		size_t num; Primitive* prim = (Primitive*)cur.leaf(num);
		PrimitiveIntersector8::intersect(ray, k, prim, num, bvh->geometry);
		rayFar = ray.tfar[k];
	      }
	  }

	  static __forceinline bool occluded1(const BVH4* bvh, NodeRef root, const size_t k, Ray8& ray,const avx3f &ray_org, const avx3f &ray_dir, const avx3f &ray_rdir, const avxf &ray_tnear, const avxf &ray_tfar, const avx3i& nearXYZ)
	  {
	    /*! stack state */
	    NodeRef stack[stackSizeSingle];  //!< stack of nodes that still need to get traversed
	    NodeRef* stackPtr = stack+1;        //!< current stack pointer
	    NodeRef* stackEnd = stack+stackSizeSingle;
	    stack[0]  = root;
      
	    /*! offsets to select the side that becomes the lower or upper bound */
	    const size_t nearX = nearXYZ.x[k];
	    const size_t nearY = nearXYZ.y[k];
	    const size_t nearZ = nearXYZ.z[k];
      
	    /*! load the ray into SIMD registers */
	    const sse3f org (ray_org .x[k],ray_org .y[k],ray_org .z[k]);
	    const sse3f rdir(ray_rdir.x[k],ray_rdir.y[k],ray_rdir.z[k]);
	    const sse3f norg = -org, org_rdir(org*rdir);
	    const ssef rayNear(ray_tnear[k]), rayFar(ray_tfar[k]); 
      
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
		    const ssef tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,rayNear));
		    const ssef tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,rayFar ));
		    const sseb vmask = cast(tNear) > cast(tFar);
		    size_t mask = movemask(vmask)^0xf;
#else
		    const ssef tNear = max(tNearX,tNearY,tNearZ,rayNear);
		    const ssef tFar  = min(tFarX ,tFarY ,tFarZ ,rayFar);
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
		    NodeRef c0 = node->child(r); unsigned int d0 = ((unsigned int*)&tNear)[r];
		    r = __bscf(mask);
		    NodeRef c1 = node->child(r); unsigned int d1 = ((unsigned int*)&tNear)[r];

#if defined(__AVX2__) && defined(ENABLE_PREFETCHING)
		    prefetchL1(((char*)c0.node()) + 0);
		    prefetchL1(((char*)c0.node()) + 64);
		    prefetchL1(((char*)c1.node()) + 0);
		    prefetchL1(((char*)c1.node()) + 64);
#endif

		    assert(c0 != BVH4::emptyNode);
		    assert(c1 != BVH4::emptyNode);
		    if (likely(mask == 0)) {
		      assert(stackPtr < stackEnd);
		      if (d0 < d1) { *stackPtr = c1; stackPtr++; cur = c0; continue; }
		      else         { *stackPtr = c0; stackPtr++; cur = c1; continue; }
		    }
		    assert(stackPtr < stackEnd);
		    stackPtr[0] = c0; 
		    assert(stackPtr < stackEnd);
		    stackPtr[1] = c1; 

		    stackPtr+=2;

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
		if (PrimitiveIntersector8::occluded(ray,k,prim,num,bvh->geometry)) {
		  ray.geomID[k] = 0;
		  return true;
		}
	      }
	    return false;
	  }


      static void intersect(avxb* valid, BVH4* bvh, Ray8& ray);
      static void occluded (avxb* valid, BVH4* bvh, Ray8& ray);
    };
  }
}

#endif
  
