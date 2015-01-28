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

#include "bvh4_intersector4_chunk.h"

#include "geometry/triangle1_intersector4_moeller.h"
#include "geometry/triangle4_intersector4_moeller.h"
#if defined (__AVX__)
#include "geometry/triangle8_intersector4_moeller.h"
#endif
#include "geometry/triangle1v_intersector4_pluecker.h"
#include "geometry/triangle4v_intersector4_pluecker.h"
#include "geometry/triangle4i_intersector4.h"
#include "geometry/virtual_accel_intersector4.h"

namespace embree
{
  namespace isa
  {
    template<typename PrimitiveIntersector4>
    void BVH4Intersector4Chunk<PrimitiveIntersector4>::intersect(sseb* valid_i, BVH4* bvh, Ray4& ray)
    {
      /* load ray */
      const sseb valid0 = *valid_i;
      const sse3f rdir = rcp_safe(ray.dir);
      const sse3f org(ray.org), org_rdir = org * rdir;
      ssef ray_tnear = select(valid0,ray.tnear,ssef(pos_inf));
      ssef ray_tfar  = select(valid0,ray.tfar ,ssef(neg_inf));
      const ssef inf = ssef(pos_inf);
      
      /* allocate stack and push root node */
      ssef    stack_near[stackSize];
      NodeRef stack_node[stackSize];
      stack_node[0] = BVH4::invalidNode;
      stack_near[0] = inf;
      stack_node[1] = bvh->root;
      stack_near[1] = ray_tnear; 
      NodeRef* stackEnd = stack_node+stackSize;
      NodeRef* __restrict__ sptr_node = stack_node + 2;
      ssef*    __restrict__ sptr_near = stack_near + 2;
      
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
        ssef curDist = *sptr_near;
        if (unlikely(none(ray_tfar > curDist))) 
          continue;
        
        while (1)
        {
          /* test if this is a leaf node */
          if (unlikely(curNode.isLeaf()))
            break;
          
          const sseb valid_node = ray_tfar > curDist;
          STAT3(normal.trav_nodes,1,popcnt(valid_node),4);
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
            const ssef lclipMinX = msub(node->lower_x[i],rdir.x,org_rdir.x);
            const ssef lclipMinY = msub(node->lower_y[i],rdir.y,org_rdir.y);
            const ssef lclipMinZ = msub(node->lower_z[i],rdir.z,org_rdir.z);
            const ssef lclipMaxX = msub(node->upper_x[i],rdir.x,org_rdir.x);
            const ssef lclipMaxY = msub(node->upper_y[i],rdir.y,org_rdir.y);
            const ssef lclipMaxZ = msub(node->upper_z[i],rdir.z,org_rdir.z);
#else
            const ssef lclipMinX = (node->lower_x[i] - org.x) * rdir.x;
            const ssef lclipMinY = (node->lower_y[i] - org.y) * rdir.y;
            const ssef lclipMinZ = (node->lower_z[i] - org.z) * rdir.z;
            const ssef lclipMaxX = (node->upper_x[i] - org.x) * rdir.x;
            const ssef lclipMaxY = (node->upper_y[i] - org.y) * rdir.y;
            const ssef lclipMaxZ = (node->upper_z[i] - org.z) * rdir.z;
#endif

#if defined(__SSE4_1__)
            const ssef lnearP = maxi(maxi(mini(lclipMinX, lclipMaxX), mini(lclipMinY, lclipMaxY)), mini(lclipMinZ, lclipMaxZ));
            const ssef lfarP  = mini(mini(maxi(lclipMinX, lclipMaxX), maxi(lclipMinY, lclipMaxY)), maxi(lclipMinZ, lclipMaxZ));
            const sseb lhit   = maxi(lnearP,ray_tnear) <= mini(lfarP,ray_tfar);      
#else
            const ssef lnearP = max(max(min(lclipMinX, lclipMaxX), min(lclipMinY, lclipMaxY)), min(lclipMinZ, lclipMaxZ));
            const ssef lfarP  = min(min(max(lclipMinX, lclipMaxX), max(lclipMinY, lclipMaxY)), max(lclipMinZ, lclipMaxZ));
            const sseb lhit   = max(lnearP,ray_tnear) <= min(lfarP,ray_tfar);      
#endif
            
            /* if we hit the child we choose to continue with that child if it 
               is closer than the current next child, or we push it onto the stack */
            if (likely(any(lhit)))
            {
              assert(sptr_node < stackEnd);
              const ssef childDist = select(lhit,lnearP,inf);
              const NodeRef child = node->children[i];
              assert(child != BVH4::emptyNode);
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
        const sseb valid_leaf = ray_tfar > curDist;
        STAT3(normal.trav_leaves,1,popcnt(valid_leaf),4);
        size_t items; const Primitive* prim = (Primitive*) curNode.leaf(items);
        PrimitiveIntersector4::intersect(valid_leaf,ray,prim,items,bvh->geometry);
        ray_tfar = select(valid_leaf,ray.tfar,ray_tfar);
      }
      AVX_ZERO_UPPER();
    }
    
    template<typename PrimitiveIntersector4>
    void BVH4Intersector4Chunk<PrimitiveIntersector4>::occluded(sseb* valid_i, BVH4* bvh, Ray4& ray)
    {
      /* load ray */
      const sseb valid = *valid_i;
      sseb terminated = !valid;
      const sse3f rdir = rcp_safe(ray.dir);
      const sse3f org(ray.org), org_rdir = org * rdir;
      ssef ray_tnear = select(valid,ray.tnear,ssef(pos_inf));
      ssef ray_tfar  = select(valid,ray.tfar ,ssef(neg_inf));
      const ssef inf = ssef(pos_inf);
      
      /* allocate stack and push root node */
      ssef    stack_near[stackSize];
      NodeRef stack_node[stackSize];
      stack_node[0] = BVH4::invalidNode;
      stack_near[0] = inf;
      stack_node[1] = bvh->root;
      stack_near[1] = ray_tnear; 
      NodeRef* stackEnd = stack_node+stackSize;
      NodeRef* __restrict__ sptr_node = stack_node + 2;
      ssef*    __restrict__ sptr_near = stack_near + 2;
      
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
        ssef curDist = *sptr_near;
        if (unlikely(none(ray_tfar > curDist))) 
          continue;
        
        while (1)
        {
          /* test if this is a leaf node */
          if (unlikely(curNode.isLeaf()))
            break;
          
          const sseb valid_node = ray_tfar > curDist;
          STAT3(shadow.trav_nodes,1,popcnt(valid_node),4);
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
            const ssef lclipMinX = msub(node->lower_x[i],rdir.x,org_rdir.x);
            const ssef lclipMinY = msub(node->lower_y[i],rdir.y,org_rdir.y);
            const ssef lclipMinZ = msub(node->lower_z[i],rdir.z,org_rdir.z);
            const ssef lclipMaxX = msub(node->upper_x[i],rdir.x,org_rdir.x);
            const ssef lclipMaxY = msub(node->upper_y[i],rdir.y,org_rdir.y);
            const ssef lclipMaxZ = msub(node->upper_z[i],rdir.z,org_rdir.z);
#else
            const ssef lclipMinX = (node->lower_x[i] - org.x) * rdir.x;
            const ssef lclipMinY = (node->lower_y[i] - org.y) * rdir.y;
            const ssef lclipMinZ = (node->lower_z[i] - org.z) * rdir.z;
            const ssef lclipMaxX = (node->upper_x[i] - org.x) * rdir.x;
            const ssef lclipMaxY = (node->upper_y[i] - org.y) * rdir.y;
            const ssef lclipMaxZ = (node->upper_z[i] - org.z) * rdir.z;
#endif

#if defined(__SSE4_1__)
            const ssef lnearP = maxi(maxi(mini(lclipMinX, lclipMaxX), mini(lclipMinY, lclipMaxY)), mini(lclipMinZ, lclipMaxZ));
            const ssef lfarP  = mini(mini(maxi(lclipMinX, lclipMaxX), maxi(lclipMinY, lclipMaxY)), maxi(lclipMinZ, lclipMaxZ));
            const sseb lhit   = maxi(lnearP,ray_tnear) <= mini(lfarP,ray_tfar);      
#else
            const ssef lnearP = max(max(min(lclipMinX, lclipMaxX), min(lclipMinY, lclipMaxY)), min(lclipMinZ, lclipMaxZ));
            const ssef lfarP  = min(min(max(lclipMinX, lclipMaxX), max(lclipMinY, lclipMaxY)), max(lclipMinZ, lclipMaxZ));
            const sseb lhit   = max(lnearP,ray_tnear) <= min(lfarP,ray_tfar);      
#endif
            
            /* if we hit the child we choose to continue with that child if it 
               is closer than the current next child, or we push it onto the stack */
            if (likely(any(lhit)))
            {
              assert(sptr_node < stackEnd);
              assert(child != BVH4::emptyNode);
              const ssef childDist = select(lhit,lnearP,inf);
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
        const sseb valid_leaf = ray_tfar > curDist;
        STAT3(shadow.trav_leaves,1,popcnt(valid_leaf),4);
        size_t items; const Primitive* prim = (Primitive*) curNode.leaf(items);
        terminated |= PrimitiveIntersector4::occluded(!terminated,ray,prim,items,bvh->geometry);
        if (all(terminated)) break;
        ray_tfar = select(terminated,ssef(neg_inf),ray_tfar);
      }
      store4i(valid & terminated,&ray.geomID,0);
      AVX_ZERO_UPPER();
    }
    
    DEFINE_INTERSECTOR4(BVH4Triangle1Intersector4ChunkMoeller, BVH4Intersector4Chunk<Triangle1Intersector4MoellerTrumbore>);
    DEFINE_INTERSECTOR4(BVH4Triangle4Intersector4ChunkMoeller, BVH4Intersector4Chunk<Triangle4Intersector4MoellerTrumbore>);
#if defined (__AVX__)
    DEFINE_INTERSECTOR4(BVH4Triangle8Intersector4ChunkMoeller, BVH4Intersector4Chunk<Triangle8Intersector4MoellerTrumbore>);
#endif
    DEFINE_INTERSECTOR4(BVH4Triangle1vIntersector4ChunkPluecker, BVH4Intersector4Chunk<Triangle1vIntersector4Pluecker>);
    DEFINE_INTERSECTOR4(BVH4Triangle4vIntersector4ChunkPluecker, BVH4Intersector4Chunk<Triangle4vIntersector4Pluecker>);
    DEFINE_INTERSECTOR4(BVH4Triangle4iIntersector4ChunkPluecker, BVH4Intersector4Chunk<Triangle4iIntersector4Pluecker>);
    DEFINE_INTERSECTOR4(BVH4VirtualIntersector4Chunk, BVH4Intersector4Chunk<VirtualAccelIntersector4>);
  }
}
