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

#pragma once

#include "bvh4.h"
#include "common/ray8.h"
#include "common/stack_item.h"

namespace embree
{
  namespace isa 
  {
    /*! Converts single ray traversal into packet traversal. */
    template<typename Intersector1>
    class BVH4Intersector8FromIntersector1
    {
    public:
      static void intersect(avxb* valid, BVH4* bvh, Ray8& ray);
      static void occluded (avxb* valid, BVH4* bvh, Ray8& ray);
    };

    /*! Single ray traversal for packets. */
    template<int types, bool robust, typename PrimitiveIntersector8>
    class BVH4Intersector8Single 
    {
      /* shortcuts for frequently used types */
      typedef typename PrimitiveIntersector8::Precalculations Precalculations;
      typedef typename PrimitiveIntersector8::Primitive Primitive;
      typedef typename BVH4::NodeRef NodeRef;
      typedef typename BVH4::Node Node;
      typedef StackItemT<NodeRef> StackItem;
      static const size_t stackSizeSingle = 1+3*BVH4::maxDepth;
      static const size_t stackSizeChunk = 4*BVH4::maxDepth+1;
      
    public:

      static __forceinline void intersect1(const BVH4* bvh, NodeRef root, const size_t k, Precalculations& pre, 
					   Ray8& ray, const avx3f &ray_org, const avx3f &ray_dir, const avx3f &ray_rdir, const avxf &ray_tnear, const avxf &ray_tfar, 
					   const avx3i& nearXYZ)
    {
      /*! stack state */
      StackItemInt32<NodeRef> stack[stackSizeSingle];  //!< stack of nodes 
      StackItemInt32<NodeRef>* stackPtr = stack + 1;        //!< current stack pointer
      StackItemInt32<NodeRef>* stackEnd = stack + stackSizeSingle;
      stack[0].ptr  = root;
      stack[0].dist = neg_inf;
      
      /*! load the ray into SIMD registers */
      const sse3f org(ray_org.x[k], ray_org.y[k], ray_org.z[k]);
      const sse3f dir(ray_dir.x[k], ray_dir.y[k], ray_dir.z[k]);
      const sse3f rdir(ray_rdir.x[k], ray_rdir.y[k], ray_rdir.z[k]);
      const sse3f org_rdir(org*rdir);
      ssef ray_near(ray_tnear[k]), ray_far(ray_tfar[k]);
      
      /*! offsets to select the side that becomes the lower or upper bound */
      const size_t nearX = nearXYZ.x[k];
      const size_t nearY = nearXYZ.y[k];
      const size_t nearZ = nearXYZ.z[k];

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
	  size_t mask; 
	  ssef tNear;

	  /*! stop if we found a leaf node */
	  if (unlikely(cur.isLeaf(types))) break;
	  STAT3(normal.trav_nodes,1,1,1);

	  /* process standard nodes */
          if (likely(cur.isNode(types)))
	    mask = cur.node()->intersect<robust>(nearX,nearY,nearZ,org,rdir,org_rdir,ray_near,ray_far,tNear); 

	  /* process motion blur nodes */
	  else if (likely(cur.isNodeMB(types)))
	    mask = cur.nodeMB()->intersect(nearX,nearY,nearZ,org,rdir,org_rdir,ray_near,ray_far,ray.time[k],tNear); 

	  /*! process nodes with unaligned bounds */
          else if (unlikely(cur.isUnalignedNode(types)))
            mask = cur.unalignedNode()->intersect(org,dir,ray_near,ray_far,tNear);

          /*! process nodes with unaligned bounds and motion blur */
          else if (unlikely(cur.isUnalignedNodeMB(types)))
            mask = cur.unalignedNodeMB()->intersect(org,dir,ray_near,ray_far,ray.time[k],tNear);

          /*! if no child is hit, pop next node */
	  const BVH4::BaseNode* node = cur.baseNode(types);
          if (unlikely(mask == 0))
            goto pop;
          
          /*! one child is hit, continue with that child */
	  size_t r = __bscf(mask);
	  if (likely(mask == 0)) {
            cur = node->child(r); cur.prefetch(types);
            assert(cur != BVH4::emptyNode);
            continue;
          }
          
          /*! two children are hit, push far child, and continue with closer child */
          NodeRef c0 = node->child(r); c0.prefetch(types); const unsigned int d0 = ((unsigned int*)&tNear)[r];
          r = __bscf(mask);
          NodeRef c1 = node->child(r); c1.prefetch(types); const unsigned int d1 = ((unsigned int*)&tNear)[r];
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
          NodeRef c = node->child(r); c.prefetch(types); unsigned int d = ((unsigned int*)&tNear)[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
          assert(c != BVH4::emptyNode);
          if (likely(mask == 0)) {
            sort(stackPtr[-1],stackPtr[-2],stackPtr[-3]);
            cur = (NodeRef) stackPtr[-1].ptr; stackPtr--;
            continue;
          }
          
          /*! four children are hit, push all onto stack and sort 4 stack items, continue with closest child */
          assert(stackPtr < stackEnd); 
          r = __bscf(mask);
          c = node->child(r); c.prefetch(types); d = *(unsigned int*)&tNear[r]; stackPtr->ptr = c; stackPtr->dist = d; stackPtr++;
          assert(c != BVH4::emptyNode);
          sort(stackPtr[-1],stackPtr[-2],stackPtr[-3],stackPtr[-4]);
          cur = (NodeRef) stackPtr[-1].ptr; stackPtr--;
        }
        
        /*! this is a leaf node */
	assert(cur != BVH4::emptyNode);
	STAT3(normal.trav_leaves, 1, 1, 1);
	size_t num; Primitive* prim = (Primitive*)cur.leaf(num);
	PrimitiveIntersector8::intersect(pre, ray, k, prim, num, bvh->scene);
        ray_far = ray.tfar[k];
      }
    }
    
      static __forceinline bool occluded1(const BVH4* bvh, NodeRef root, const size_t k, Precalculations& pre, 
					  Ray8& ray,const avx3f &ray_org, const avx3f &ray_dir, const avx3f &ray_rdir, const avxf &ray_tnear, const avxf &ray_tfar, 
					  const avx3i& nearXYZ)
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
      const sse3f dir(ray_dir.x[k], ray_dir.y[k], ray_dir.z[k]);
      const sse3f rdir(ray_rdir.x[k],ray_rdir.y[k],ray_rdir.z[k]);
      const sse3f norg = -org, org_rdir(org*rdir);
      const ssef ray_near(ray_tnear[k]), ray_far(ray_tfar[k]); 
      
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
	  size_t mask; 
	  ssef tNear;

	  /*! stop if we found a leaf node */
	  if (unlikely(cur.isLeaf(types))) break;
	  STAT3(shadow.trav_nodes,1,1,1);

	  /* process standard nodes */
          if (likely(cur.isNode(types)))
	    mask = cur.node()->intersect<robust>(nearX,nearY,nearZ,org,rdir,org_rdir,ray_near,ray_far,tNear); 

	  /* process motion blur nodes */
	  else if (likely(cur.isNodeMB(types)))
	    mask = cur.nodeMB()->intersect(nearX,nearY,nearZ,org,rdir,org_rdir,ray_near,ray_far,ray.time[k],tNear); 

	  /*! process nodes with unaligned bounds */
          else if (unlikely(cur.isUnalignedNode(types)))
            mask = cur.unalignedNode()->intersect(org,dir,ray_near,ray_far,tNear);

          /*! process nodes with unaligned bounds and motion blur */
          else if (unlikely(cur.isUnalignedNodeMB(types)))
            mask = cur.unalignedNodeMB()->intersect(org,dir,ray_near,ray_far,ray.time[k],tNear);
	  
          /*! if no child is hit, pop next node */
	  const BVH4::BaseNode* node = cur.baseNode(types);
          if (unlikely(mask == 0))
            goto pop;
	  
	  /*! one child is hit, continue with that child */
          size_t r = __bscf(mask);
          if (likely(mask == 0)) {
            cur = node->child(r); cur.prefetch(types); 
            assert(cur != BVH4::emptyNode);
            continue;
          }
          
          /*! two children are hit, push far child, and continue with closer child */
          NodeRef c0 = node->child(r); c0.prefetch(types); const unsigned int d0 = ((unsigned int*)&tNear)[r];
          r = __bscf(mask);
          NodeRef c1 = node->child(r); c1.prefetch(types); const unsigned int d1 = ((unsigned int*)&tNear)[r];
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
          cur = node->child(r); cur.prefetch(types);
          assert(cur != BVH4::emptyNode);
          if (likely(mask == 0)) continue;
          assert(stackPtr < stackEnd);
          *stackPtr = cur; stackPtr++;
          
          /*! four children are hit */
          cur = node->child(3); cur.prefetch(types);
          assert(cur != BVH4::emptyNode);
        }
	
	/*! this is a leaf node */
	assert(cur != BVH4::emptyNode);
	STAT3(shadow.trav_leaves,1,1,1);
	size_t num; Primitive* prim = (Primitive*) cur.leaf(num);
	if (PrimitiveIntersector8::occluded(pre,ray,k,prim,num,bvh->scene)) {
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
