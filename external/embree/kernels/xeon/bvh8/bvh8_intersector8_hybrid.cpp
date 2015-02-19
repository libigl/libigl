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

#include "bvh8_intersector8_hybrid.h"
#include "geometry/triangle4_intersector8_moeller.h"
#include "geometry/triangle8_intersector8_moeller.h"

#define DBG(x) 

#define SWITCH_THRESHOLD 7


namespace embree
{
  namespace isa
  { 

    template<typename PrimitiveIntersector8>
    __forceinline void BVH8Intersector8Hybrid<PrimitiveIntersector8>::intersect1(const BVH8* bvh, NodeRef root, const size_t k, Precalculations& pre, Ray8& ray,const avx3f &ray_org, const avx3f &ray_dir, const avx3f &ray_rdir, const avxf &ray_tnear, const avxf &ray_tfar, const avx3i& nearXYZ)
    {
      /*! stack state */
      StackItemInt32<NodeRef> stack[stackSizeSingle];  //!< stack of nodes 
      StackItemInt32<NodeRef>* stackPtr = stack+1;        //!< current stack pointer
      StackItemInt32<NodeRef>* stackEnd = stack+stackSizeSingle;
      stack[0].ptr = root;
      stack[0].dist = neg_inf;
      
      /*! offsets to select the side that becomes the lower or upper bound */
      const size_t nearX = nearXYZ.x[k];
      const size_t nearY = nearXYZ.y[k];
      const size_t nearZ = nearXYZ.z[k];

      /*! load the ray into SIMD registers */
      const avx3f org (ray_org .x[k],ray_org .y[k],ray_org .z[k]);
      const avx3f rdir(ray_rdir.x[k],ray_rdir.y[k],ray_rdir.z[k]);
      const avx3f norg = -org, org_rdir(org*rdir);
      avxf rayNear(ray_tnear[k]), rayFar(ray_tfar[k]);
     
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
          const avxf tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,rayNear));
          const avxf tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,rayFar ));
          const avxb vmask = cast(tNear) > cast(tFar);
          size_t mask = movemask(vmask)^0xff;
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
        size_t num; Triangle* prim = (Triangle*) cur.leaf(num);
        PrimitiveIntersector8::intersect(pre,ray,k,prim,num,bvh->scene);
        rayFar = ray.tfar[k];
      }
    }
   
    
    template<typename PrimitiveIntersector8>    
    void BVH8Intersector8Hybrid<PrimitiveIntersector8>::intersect(avxb* valid_i, BVH8* bvh, Ray8& ray)
    {
      /* load ray */
      const avxb valid0 = *valid_i;
      avx3f ray_org = ray.org;
      avx3f ray_dir = ray.dir;
      avxf ray_tnear = ray.tnear, ray_tfar  = ray.tfar;

      const avx3f rdir = rcp_safe(ray_dir);
      const avx3f org(ray_org), org_rdir = org * rdir;
      ray_tnear = select(valid0,ray_tnear,avxf(pos_inf));
      ray_tfar  = select(valid0,ray_tfar ,avxf(neg_inf));
      const avxf inf = avxf(pos_inf);
      Precalculations pre(valid0,ray);

      /* compute near/far per ray */
      avx3i nearXYZ;
      nearXYZ.x = select(rdir.x >= 0.0f,avxi(0*(int)sizeof(avxf)),avxi(1*(int)sizeof(avxf)));
      nearXYZ.y = select(rdir.y >= 0.0f,avxi(2*(int)sizeof(avxf)),avxi(3*(int)sizeof(avxf)));
      nearXYZ.z = select(rdir.z >= 0.0f,avxi(4*(int)sizeof(avxf)),avxi(5*(int)sizeof(avxf)));

      /* allocate stack and push root node */
      avxf    stack_near[stackSizeChunk];
      NodeRef stack_node[stackSizeChunk];
      stack_node[0] = BVH8::invalidNode;
      stack_near[0] = inf;
      stack_node[1] = bvh->root;
      stack_near[1] = ray_tnear; 
      NodeRef* stackEnd = stack_node+stackSizeChunk;
      NodeRef* __restrict__ sptr_node = stack_node + 2;
      avxf*    __restrict__ sptr_near = stack_near + 2;
      
      while (1)
      {
        /* pop next node from stack */
        assert(sptr_node > stack_node);
        sptr_node--;
        sptr_near--;
        NodeRef cur = *sptr_node;
        if (unlikely(cur == BVH8::invalidNode)) {
          assert(sptr_node == stack_node);
          break;
        }
        
        /* cull node if behind closest hit point */
        avxf curDist = *sptr_near;
        const avxb active = curDist < ray_tfar;
        if (unlikely(none(active)))
          continue;
        
        /* switch to single ray traversal */
#if !defined(__WIN32__) || defined(__X86_64__)
        size_t bits = movemask(active);
        if (unlikely(__popcnt(bits) <= SWITCH_THRESHOLD)) {
          for (size_t i=__bsf(bits); bits!=0; bits=__btc(bits,i), i=__bsf(bits)) {
            intersect1(bvh,cur,i,pre,ray,ray_org,ray_dir,rdir,ray_tnear,ray_tfar,nearXYZ);
          }
          ray_tfar = ray.tfar;
          continue;
        }
#endif

        while (1)
        {
          /* test if this is a leaf node */
          if (unlikely(cur.isLeaf()))
            break;
          
          const avxb valid_node = ray_tfar > curDist;
          STAT3(normal.trav_nodes,1,popcnt(valid_node),8);
          const Node* __restrict__ const node = (Node*)cur.node();
          
          /* pop of next node */
          assert(sptr_node > stack_node);
          sptr_node--;
          sptr_near--;
          cur = *sptr_node;
          curDist = *sptr_near;
          
          for (unsigned i=0; i<BVH8::N; i++)
          {
            const NodeRef child = node->children[i];
            if (unlikely(child == BVH8::emptyNode)) break;
            
#if defined(__AVX2__)
            const avxf lclipMinX = msub(node->lower_x[i],rdir.x,org_rdir.x);
            const avxf lclipMinY = msub(node->lower_y[i],rdir.y,org_rdir.y);
            const avxf lclipMinZ = msub(node->lower_z[i],rdir.z,org_rdir.z);
            const avxf lclipMaxX = msub(node->upper_x[i],rdir.x,org_rdir.x);
            const avxf lclipMaxY = msub(node->upper_y[i],rdir.y,org_rdir.y);
            const avxf lclipMaxZ = msub(node->upper_z[i],rdir.z,org_rdir.z);
            const avxf lnearP = maxi(maxi(mini(lclipMinX, lclipMaxX), mini(lclipMinY, lclipMaxY)), mini(lclipMinZ, lclipMaxZ));
            const avxf lfarP  = mini(mini(maxi(lclipMinX, lclipMaxX), maxi(lclipMinY, lclipMaxY)), maxi(lclipMinZ, lclipMaxZ));
            const avxb lhit   = maxi(lnearP,ray_tnear) <= mini(lfarP,ray_tfar);      
#else
            const avxf lclipMinX = (node->lower_x[i] - org.x) * rdir.x;
            const avxf lclipMinY = (node->lower_y[i] - org.y) * rdir.y;
            const avxf lclipMinZ = (node->lower_z[i] - org.z) * rdir.z;
            const avxf lclipMaxX = (node->upper_x[i] - org.x) * rdir.x;
            const avxf lclipMaxY = (node->upper_y[i] - org.y) * rdir.y;
            const avxf lclipMaxZ = (node->upper_z[i] - org.z) * rdir.z;
            const avxf lnearP = max(max(min(lclipMinX, lclipMaxX), min(lclipMinY, lclipMaxY)), min(lclipMinZ, lclipMaxZ));
            const avxf lfarP  = min(min(max(lclipMinX, lclipMaxX), max(lclipMinY, lclipMaxY)), max(lclipMinZ, lclipMaxZ));
            const avxb lhit   = max(lnearP,ray_tnear) <= min(lfarP,ray_tfar);      
#endif
            
            /* if we hit the child we choose to continue with that child if it 
               is closer than the current next child, or we push it onto the stack */
            if (likely(any(lhit)))
            {
              assert(sptr_node < stackEnd);
              const avxf childDist = select(lhit,lnearP,inf);
              const NodeRef child = node->children[i];
              assert(child != BVH8::emptyNode);
              
              /* push cur node onto stack and continue with hit child */
              if (any(childDist < curDist))
              {
                *sptr_node = cur;
                *sptr_near = curDist; 
                curDist = childDist;
                cur = child;
		sptr_node++;
		sptr_near++;
              }              
              /* push hit child onto stack */
              else {
                *sptr_node = child;
                *sptr_near = childDist;               
		sptr_node++;
		sptr_near++;

              }
            }	      
          }
        }
        
        /* return if stack is empty */
        if (unlikely(cur == BVH8::invalidNode)) {
          assert(sptr_node == stack_node);
          break;
        }
        
        /* intersect leaf */
	assert(cur != BVH8::emptyNode);
        const avxb valid_leaf = ray_tfar > curDist;
        STAT3(normal.trav_leaves,1,popcnt(valid_leaf),8);
        size_t items; const Triangle* prim = (Triangle*) cur.leaf(items);
        PrimitiveIntersector8::intersect(valid_leaf,pre,ray,prim,items,bvh->scene);
        ray_tfar = select(valid_leaf,ray.tfar,ray_tfar);
      }
      AVX_ZERO_UPPER();
    }
    


    template<typename PrimitiveIntersector8>
    __forceinline bool BVH8Intersector8Hybrid<PrimitiveIntersector8>::occluded1(const BVH8* bvh, NodeRef root, const size_t k, Precalculations& pre, Ray8& ray,const avx3f &ray_org, const avx3f &ray_dir, const avx3f &ray_rdir, const avxf &ray_tnear, const avxf &ray_tfar, const avx3i& nearXYZ)
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
      const avx3f org (ray_org .x[k],ray_org .y[k],ray_org .z[k]);
      const avx3f rdir(ray_rdir.x[k],ray_rdir.y[k],ray_rdir.z[k]);
      const avx3f norg = -org, org_rdir(org*rdir);
      const avxf rayNear(ray_tnear[k]), rayFar(ray_tfar[k]); 

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
          const avxf tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,rayNear));
          const avxf tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,rayFar ));
          const avxb vmask = cast(tNear) > cast(tFar);
          size_t mask = movemask(vmask)^0xff;
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
        size_t num; Triangle* prim = (Triangle*) cur.leaf(num);
        if (PrimitiveIntersector8::occluded(pre,ray,k,prim,num,bvh->scene)) {
          //ray.geomID = 0;
          //break;
	  return true;
        }
      }
      return false;
    }

     template<typename PrimitiveIntersector8>
    void BVH8Intersector8Hybrid<PrimitiveIntersector8>::occluded(avxb* valid_i, BVH8* bvh, Ray8& ray)
    {
      /* load ray */
      const avxb valid = *valid_i;
      avxb terminated = !valid;
      avx3f ray_org = ray.org, ray_dir = ray.dir;
      avxf ray_tnear = ray.tnear, ray_tfar  = ray.tfar;
      const avx3f rdir = rcp_safe(ray_dir);
      const avx3f org(ray_org), org_rdir = org * rdir;
      ray_tnear = select(valid,ray_tnear,avxf(pos_inf));
      ray_tfar  = select(valid,ray_tfar ,avxf(neg_inf));
      const avxf inf = avxf(pos_inf);
      Precalculations pre(valid,ray);

      /* compute near/far per ray */
      avx3i nearXYZ;
      nearXYZ.x = select(rdir.x >= 0.0f,avxi(0*(int)sizeof(avxf)),avxi(1*(int)sizeof(avxf)));
      nearXYZ.y = select(rdir.y >= 0.0f,avxi(2*(int)sizeof(avxf)),avxi(3*(int)sizeof(avxf)));
      nearXYZ.z = select(rdir.z >= 0.0f,avxi(4*(int)sizeof(avxf)),avxi(5*(int)sizeof(avxf)));

      /* allocate stack and push root node */
      avxf    stack_near[stackSizeChunk];
      NodeRef stack_node[stackSizeChunk];
      stack_node[0] = BVH8::invalidNode;
      stack_near[0] = inf;
      stack_node[1] = bvh->root;
      stack_near[1] = ray_tnear; 
      NodeRef* stackEnd = stack_node+stackSizeChunk;
      NodeRef* __restrict__ sptr_node = stack_node + 2;
      avxf*    __restrict__ sptr_near = stack_near + 2;

      while (1)
      {
        /* pop next node from stack */
        assert(sptr_node > stack_node);
        sptr_node--;
        sptr_near--;
        NodeRef cur = *sptr_node;
        if (unlikely(cur == BVH8::invalidNode)) {
          assert(sptr_node == stack_node);
          break;
        }

        /* cull node if behind closest hit point */
        avxf curDist = *sptr_near;
        const avxb active = curDist < ray_tfar;
        if (unlikely(none(active))) 
          continue;
        
        /* switch to single ray traversal */
#if !defined(__WIN32__) || defined(__X86_64__)
        size_t bits = movemask(active);
        if (unlikely(__popcnt(bits) <= SWITCH_THRESHOLD)) {
          for (size_t i=__bsf(bits); bits!=0; bits=__btc(bits,i), i=__bsf(bits)) {
            if (occluded1(bvh,cur,i,pre,ray,ray_org,ray_dir,rdir,ray_tnear,ray_tfar,nearXYZ))
              terminated[i] = -1;
          }
          if (all(terminated)) break;
          ray_tfar = select(terminated,avxf(neg_inf),ray_tfar);
          continue;
        }
#endif
                
        while (1)
        {
          /* test if this is a leaf node */
          if (unlikely(cur.isLeaf()))
            break;
          
          const avxb valid_node = ray_tfar > curDist;
          STAT3(shadow.trav_nodes,1,popcnt(valid_node),8);
          const Node* __restrict__ const node = (Node*)cur.node();
          
          /* pop of next node */
          assert(sptr_node > stack_node);
          sptr_node--;
          sptr_near--;
          cur = *sptr_node;
          curDist = *sptr_near;
          
          for (unsigned i=0; i<BVH8::N; i++)
          {
            const NodeRef child = node->children[i];
            if (unlikely(child == BVH8::emptyNode)) break;
            
#if defined(__AVX2__)
            const avxf lclipMinX = msub(node->lower_x[i],rdir.x,org_rdir.x);
            const avxf lclipMinY = msub(node->lower_y[i],rdir.y,org_rdir.y);
            const avxf lclipMinZ = msub(node->lower_z[i],rdir.z,org_rdir.z);
            const avxf lclipMaxX = msub(node->upper_x[i],rdir.x,org_rdir.x);
            const avxf lclipMaxY = msub(node->upper_y[i],rdir.y,org_rdir.y);
            const avxf lclipMaxZ = msub(node->upper_z[i],rdir.z,org_rdir.z);
            const avxf lnearP = maxi(maxi(mini(lclipMinX, lclipMaxX), mini(lclipMinY, lclipMaxY)), mini(lclipMinZ, lclipMaxZ));
            const avxf lfarP  = mini(mini(maxi(lclipMinX, lclipMaxX), maxi(lclipMinY, lclipMaxY)), maxi(lclipMinZ, lclipMaxZ));
            const avxb lhit   = maxi(lnearP,ray_tnear) <= mini(lfarP,ray_tfar);      
#else
            const avxf lclipMinX = (node->lower_x[i] - org.x) * rdir.x;
            const avxf lclipMinY = (node->lower_y[i] - org.y) * rdir.y;
            const avxf lclipMinZ = (node->lower_z[i] - org.z) * rdir.z;
            const avxf lclipMaxX = (node->upper_x[i] - org.x) * rdir.x;
            const avxf lclipMaxY = (node->upper_y[i] - org.y) * rdir.y;
            const avxf lclipMaxZ = (node->upper_z[i] - org.z) * rdir.z;
            const avxf lnearP = max(max(min(lclipMinX, lclipMaxX), min(lclipMinY, lclipMaxY)), min(lclipMinZ, lclipMaxZ));
            const avxf lfarP  = min(min(max(lclipMinX, lclipMaxX), max(lclipMinY, lclipMaxY)), max(lclipMinZ, lclipMaxZ));
            const avxb lhit   = max(lnearP,ray_tnear) <= min(lfarP,ray_tfar);      
#endif
            
            /* if we hit the child we choose to continue with that child if it 
               is closer than the current next child, or we push it onto the stack */
            if (likely(any(lhit)))
            {
              assert(sptr_node < stackEnd);
              assert(child != BVH8::emptyNode);
              const avxf childDist = select(lhit,lnearP,inf);
              sptr_node++;
              sptr_near++;
              
              /* push cur node onto stack and continue with hit child */
              if (any(childDist < curDist))
              {
                *(sptr_node-1) = cur;
                *(sptr_near-1) = curDist; 
                curDist = childDist;
                cur = child;
              }
              
              /* push hit child onto stack */
              else {
                *(sptr_node-1) = child;
                *(sptr_near-1) = childDist; 
              }
            }	      
          }
        }
        
        /* return if stack is empty */
        if (unlikely(cur == BVH8::invalidNode)) {
          assert(sptr_node == stack_node);
          break;
        }
        
        /* intersect leaf */
	assert(cur != BVH8::emptyNode);
        const avxb valid_leaf = ray_tfar > curDist;
        STAT3(shadow.trav_leaves,1,popcnt(valid_leaf),8);
        size_t items; const Triangle* prim = (Triangle*) cur.leaf(items);
        terminated |= PrimitiveIntersector8::occluded(!terminated,pre,ray,prim,items,bvh->scene);
        if (all(terminated)) break;
        ray_tfar = select(terminated,avxf(neg_inf),ray_tfar);
      }
      store8i(valid & terminated,&ray.geomID,0);
      AVX_ZERO_UPPER();
    }

    DEFINE_INTERSECTOR8(BVH8Triangle4Intersector8HybridMoeller,BVH8Intersector8Hybrid<LeafIterator8_1<Triangle4Intersector8MoellerTrumbore<LeafMode COMMA true> > >);
    DEFINE_INTERSECTOR8(BVH8Triangle4Intersector8HybridMoellerNoFilter,BVH8Intersector8Hybrid<LeafIterator8_1<Triangle4Intersector8MoellerTrumbore<LeafMode COMMA false> > >);
    
    DEFINE_INTERSECTOR8(BVH8Triangle8Intersector8HybridMoeller,BVH8Intersector8Hybrid<LeafIterator8_1<Triangle8Intersector8MoellerTrumbore<LeafMode COMMA true> > >);
    DEFINE_INTERSECTOR8(BVH8Triangle8Intersector8HybridMoellerNoFilter,BVH8Intersector8Hybrid<LeafIterator8_1<Triangle8Intersector8MoellerTrumbore<LeafMode COMMA false> > >);
  }
}  

