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

#include "bvh4_intersector8_hybrid.h"

#include "geometry/triangle4_intersector8_moeller.h"
#include "geometry/triangle8_intersector8_moeller.h"
#include "geometry/triangle4v_intersector8_pluecker.h"

#define SWITCH_THRESHOLD 6


namespace embree
{
  namespace isa
  {
    template<typename PrimitiveIntersector8>
    void BVH4Intersector8Hybrid<PrimitiveIntersector8>::intersect(avxb* valid_i, BVH4* bvh, Ray8& ray)
    {
      /* load ray */
      const avxb valid0 = *valid_i;
      avx3f ray_org = ray.org;
      avx3f ray_dir = ray.dir;
      avxf ray_tnear = ray.tnear, ray_tfar  = ray.tfar;
#if defined(__FIX_RAYS__)
      const avxf float_range = 1.8E19;
      ray_org = clamp(ray_org,avx3f(-float_range),avx3f(+float_range));
      ray_dir = clamp(ray_dir,avx3f(-float_range),avx3f(+float_range));
      ray_tnear = max(ray_tnear,FLT_MIN); 
      ray_tfar  = min(ray_tfar,float(inf)); 
#endif
      const avx3f rdir = rcp_safe(ray_dir);
      const avx3f org(ray_org), org_rdir = org * rdir;
      ray_tnear = select(valid0,ray_tnear,avxf(pos_inf));
      ray_tfar  = select(valid0,ray_tfar ,avxf(neg_inf));
      const avxf inf = avxf(pos_inf);

      /* compute near/far per ray */
      avx3i nearXYZ;
      nearXYZ.x = select(rdir.x >= 0.0f,avxi(0*(int)sizeof(ssef)),avxi(1*(int)sizeof(ssef)));
      nearXYZ.y = select(rdir.y >= 0.0f,avxi(2*(int)sizeof(ssef)),avxi(3*(int)sizeof(ssef)));
      nearXYZ.z = select(rdir.z >= 0.0f,avxi(4*(int)sizeof(ssef)),avxi(5*(int)sizeof(ssef)));

      /* allocate stack and push root node */
      avxf    stack_near[stackSizeChunk];
      NodeRef stack_node[stackSizeChunk];
      stack_node[0] = BVH4::invalidNode;
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
        NodeRef curNode = *sptr_node;
        if (unlikely(curNode == BVH4::invalidNode)) {
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
            intersect1(bvh, curNode, i, ray, ray_org, ray_dir, rdir, ray_tnear, ray_tfar, nearXYZ);
          }
          ray_tfar = ray.tfar;
          continue;
        }
#endif

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
          
          for (size_t i=0; i<4; i++)
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
              const avxf childDist = select(lhit,lnearP,inf);
              const NodeRef child = node->children[i];
              assert(child != BVH4::emptyNode);

#if defined(__AVX2__) && defined(ENABLE_PREFETCHING)
	      prefetchL1(((char*)child.node()) + 0);
	      prefetchL1(((char*)child.node()) + 64);
#endif              
              
              /* push cur node onto stack and continue with hit child */
              if (any(childDist < curDist))
              {
                *sptr_node = curNode;
                *sptr_near = curDist; 
                curDist = childDist;
                curNode = child;
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
    void BVH4Intersector8Hybrid<PrimitiveIntersector8>::occluded(avxb* valid_i, BVH4* bvh, Ray8& ray)
    {
      /* load ray */
      const avxb valid = *valid_i;
      avxb terminated = !valid;
      avx3f ray_org = ray.org, ray_dir = ray.dir;
      avxf ray_tnear = ray.tnear, ray_tfar  = ray.tfar;
#if defined(__FIX_RAYS__)
      const avxf float_range = 0.1f*FLT_MAX;
      ray_org = clamp(ray_org,avx3f(-float_range),avx3f(+float_range));
      ray_dir = clamp(ray_dir,avx3f(-float_range),avx3f(+float_range));
      ray_tnear = max(ray_tnear,FLT_MIN); 
      ray_tfar  = min(ray_tfar,float(inf)); 
#endif
      const avx3f rdir = rcp_safe(ray_dir);
      const avx3f org(ray_org), org_rdir = org * rdir;
      ray_tnear = select(valid,ray_tnear,avxf(pos_inf));
      ray_tfar  = select(valid,ray_tfar ,avxf(neg_inf));
      const avxf inf = avxf(pos_inf);
      
      /* compute near/far per ray */
      avx3i nearXYZ;
      nearXYZ.x = select(rdir.x >= 0.0f,avxi(0*(int)sizeof(ssef)),avxi(1*(int)sizeof(ssef)));
      nearXYZ.y = select(rdir.y >= 0.0f,avxi(2*(int)sizeof(ssef)),avxi(3*(int)sizeof(ssef)));
      nearXYZ.z = select(rdir.z >= 0.0f,avxi(4*(int)sizeof(ssef)),avxi(5*(int)sizeof(ssef)));

      /* allocate stack and push root node */
      avxf    stack_near[stackSizeChunk];
      NodeRef stack_node[stackSizeChunk];
      stack_node[0] = BVH4::invalidNode;
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
        NodeRef curNode = *sptr_node;
        if (unlikely(curNode == BVH4::invalidNode)) {
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
            if (occluded1(bvh,curNode,i,ray,ray_org,ray_dir,rdir,ray_tnear,ray_tfar,nearXYZ))
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
          
          for (size_t i=0; i<4; i++)
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

#if defined(__AVX2__) && defined(ENABLE_PREFETCHING)
	      prefetchL1(((char*)child.node()) + 0);
	      prefetchL1(((char*)child.node()) + 64);
#endif              
              
              /* push cur node onto stack and continue with hit child */
              if (any(childDist < curDist))
		{
		  *sptr_node = curNode;
		  *sptr_near = curDist; 
		  curDist = childDist;
		  curNode = child;
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
        if (unlikely(curNode == BVH4::invalidNode)) {
          assert(sptr_node == stack_node);
          break;
        }
        
        /* intersect leaf */
        const avxb valid_leaf = ray_tfar > curDist;
        STAT3(shadow.trav_leaves,1,popcnt(valid_leaf),8);
        size_t items; const Primitive* prim = (Primitive*) curNode.leaf(items);
        terminated |= PrimitiveIntersector8::occluded(!terminated,ray,prim,items,bvh->geometry);
        if (all(terminated)) break;
        ray_tfar = select(terminated,avxf(neg_inf),ray_tfar);
      }
      store8i(valid & terminated,&ray.geomID,0);
      AVX_ZERO_UPPER();
    }
    
    DEFINE_INTERSECTOR8(BVH4Triangle4Intersector8HybridMoeller, BVH4Intersector8Hybrid<Triangle4Intersector8MoellerTrumbore>);
    DEFINE_INTERSECTOR8(BVH4Triangle8Intersector8HybridMoeller, BVH4Intersector8Hybrid<Triangle8Intersector8MoellerTrumbore>);
    DEFINE_INTERSECTOR8(BVH4Triangle4vIntersector8HybridPluecker, BVH4Intersector8Hybrid<Triangle4vIntersector8Pluecker>);

  }
}
