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

#include "bvh4_intersector8_chunk.h"

#include "geometry/triangle1_intersector8_moeller.h"
#include "geometry/triangle4_intersector8_moeller.h"
#include "geometry/triangle8_intersector8_moeller.h"
#include "geometry/triangle1v_intersector8_pluecker.h"
#include "geometry/triangle4v_intersector8_pluecker.h"
#include "geometry/triangle4i_intersector8.h"
#include "geometry/virtual_accel_intersector8.h"

namespace embree
{
  namespace isa
  {
    template<typename PrimitiveIntersector8>
    void BVH4Intersector8Chunk<PrimitiveIntersector8>::intersect(avxb* valid_i, BVH4* bvh, Ray8& ray)
    {
      /* load ray */
      const avxb valid0 = *valid_i;
      const avx3f rdir = rcp_safe(ray.dir);
      const avx3f org(ray.org), org_rdir = org * rdir;
      avxf ray_tnear = select(valid0,ray.tnear,pos_inf);
      avxf ray_tfar  = select(valid0,ray.tfar ,neg_inf);
      const avxf inf = avxf(pos_inf);
      
      /* allocate stack and push root node */
      avxf    stack_near[stackSize];
      NodeRef stack_node[stackSize];
      stack_node[0] = BVH4::invalidNode;
      stack_near[0] = inf;
      stack_node[1] = bvh->root;
      stack_near[1] = ray_tnear; 
      NodeRef* stackEnd = stack_node+stackSize;
      NodeRef* __restrict__ sptr_node = stack_node + 2;
      avxf*    __restrict__ sptr_near = stack_near + 2;
      
      while (1)
      {
        /* pop next node from stack */
        assert(sptr_node > stack_node);
        sptr_node--;
        sptr_near--;
        NodeRef curNode = *sptr_node;
        if (unlikely(curNode == BVH4::invalidNode)) {
          assert(sptr_node == stack_node);
          break;
        }
        
        /* cull node if behind closest hit point */
        avxf curDist = *sptr_near;
        if (unlikely(none(ray_tfar > curDist))) 
          continue;
        
        while (1)
        {
          /* test if this is a leaf node */
          if (unlikely(curNode.isLeaf()))
            break;
          
          const avxb valid_node = ray_tfar > curDist;
          STAT3(normal.trav_nodes,1,popcnt(valid_node),8);
          const Node* __restrict__ const node = curNode.node();
          
          /* pop of next node */
          assert(sptr_node > stack_node);
          sptr_node--;
          sptr_near--;
          curNode = *sptr_node; 
          curDist = *sptr_near;
          
#pragma unroll(4)
          for (unsigned i=0; i<4; i++)
          {
            const NodeRef child = node->children[i];
            if (unlikely(child == BVH4::emptyNode)) break;
            
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
              assert(child != BVH4::emptyNode);
              const avxf childDist = select(lhit,lnearP,inf);
              sptr_node++;
              sptr_near++;
              
              /* push cur node onto stack and continue with hit child */
              if (any(childDist < curDist))
              {
                *(sptr_node-1) = curNode;
                *(sptr_near-1) = curDist; 
                curDist = childDist;
                curNode = child;
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
        if (unlikely(curNode == BVH4::invalidNode)) {
          assert(sptr_node == stack_node);
          break;
        }
        
        /* intersect leaf */
        const avxb valid_leaf = ray_tfar > curDist;
        STAT3(normal.trav_leaves,1,popcnt(valid_leaf),8);
        size_t items; const Primitive* prim = (Primitive*) curNode.leaf(items);
        PrimitiveIntersector8::intersect(valid_leaf,ray,prim,items,bvh->geometry);
        ray_tfar = select(valid_leaf,ray.tfar,ray_tfar);
      }
      AVX_ZERO_UPPER();
    }
    
    template<typename PrimitiveIntersector8>
    void BVH4Intersector8Chunk<PrimitiveIntersector8>::occluded(avxb* valid_i, BVH4* bvh, Ray8& ray)
    {
      /* load ray */
      const avxb valid = *valid_i;
      avxb terminated = !valid;
      const avx3f rdir = rcp_safe(ray.dir);
      const avx3f org(ray.org), org_rdir = org * rdir;
      avxf ray_tnear = select(valid,ray.tnear,pos_inf);
      avxf ray_tfar  = select(valid,ray.tfar ,neg_inf);
      const avxf inf = avxf(pos_inf);
      
      /* allocate stack and push root node */
      avxf    stack_near[stackSize];
      NodeRef stack_node[stackSize];
      stack_node[0] = BVH4::invalidNode;
      stack_near[0] = inf;
      stack_node[1] = bvh->root;
      stack_near[1] = ray_tnear; 
      NodeRef* stackEnd = stack_node+stackSize;
      NodeRef* __restrict__ sptr_node = stack_node + 2;
      avxf*    __restrict__ sptr_near = stack_near + 2;
      
      while (1)
      {
        /* pop next node from stack */
        assert(sptr_node > stack_node);
        sptr_node--;
        sptr_near--;
        NodeRef curNode = *sptr_node;
        if (unlikely(curNode == BVH4::invalidNode)) {
          assert(sptr_node == stack_node);
          break;
        }
        
        /* cull node if behind closest hit point */
        avxf curDist = *sptr_near;
        if (unlikely(none(ray_tfar > curDist))) 
          continue;
        
        while (1)
        {
          /* test if this is a leaf node */
          if (unlikely(curNode.isLeaf()))
            break;
          
          const avxb valid_node = ray_tfar > curDist;
          STAT3(shadow.trav_nodes,1,popcnt(valid_node),8);
          const Node* __restrict__ const node = curNode.node();
          
          /* pop of next node */
          assert(sptr_node > stack_node);
          sptr_node--;
          sptr_near--;
          curNode = *sptr_node;
          curDist = *sptr_near;
          
#pragma unroll(4)
          for (unsigned i=0; i<4; i++)
          {
            const NodeRef child = node->children[i];
            if (unlikely(child == BVH4::emptyNode)) break;
            
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
              assert(child != BVH4::emptyNode);
              const avxf childDist = select(lhit,lnearP,inf);
              sptr_node++;
              sptr_near++;
              
              /* push cur node onto stack and continue with hit child */
              if (any(childDist < curDist))
              {
                *(sptr_node-1) = curNode;
                *(sptr_near-1) = curDist; 
                curDist = childDist;
                curNode = child;
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
        if (unlikely(curNode == BVH4::invalidNode)) {
          assert(sptr_node == stack_node);
          break;
        }
        
        /* intersect leaf */
        const avxb valid_leaf = ray_tfar > curDist;
        STAT3(shadow.trav_leaves,1,popcnt(valid_leaf),8);
        size_t items; const Primitive* prim = (Primitive*) curNode.leaf(items);
        terminated |= valid_leaf & PrimitiveIntersector8::occluded(valid_leaf,ray,prim,items,bvh->geometry);
        if (all(terminated)) break;
        ray_tfar = select(terminated,neg_inf,ray_tfar);
      }
      store8i(valid & terminated,&ray.geomID,0);
      AVX_ZERO_UPPER();
    }
    
    DEFINE_INTERSECTOR8(BVH4Triangle1Intersector8ChunkMoeller, BVH4Intersector8Chunk<Triangle1Intersector8MoellerTrumbore>);
    DEFINE_INTERSECTOR8(BVH4Triangle4Intersector8ChunkMoeller, BVH4Intersector8Chunk<Triangle4Intersector8MoellerTrumbore>);
    DEFINE_INTERSECTOR8(BVH4Triangle8Intersector8ChunkMoeller, BVH4Intersector8Chunk<Triangle8Intersector8MoellerTrumbore>);
    DEFINE_INTERSECTOR8(BVH4Triangle1vIntersector8ChunkPluecker, BVH4Intersector8Chunk<Triangle1vIntersector8Pluecker>);
    DEFINE_INTERSECTOR8(BVH4Triangle4vIntersector8ChunkPluecker, BVH4Intersector8Chunk<Triangle4vIntersector8Pluecker>);
    DEFINE_INTERSECTOR8(BVH4Triangle4iIntersector8ChunkPluecker, BVH4Intersector8Chunk<Triangle4iIntersector8Pluecker>);
    DEFINE_INTERSECTOR8(BVH4VirtualIntersector8Chunk, BVH4Intersector8Chunk<VirtualAccelIntersector8>);
  }
}
