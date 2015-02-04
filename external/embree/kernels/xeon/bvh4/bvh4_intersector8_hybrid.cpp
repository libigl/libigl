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

#include "bvh4_intersector8_hybrid.h"
#include "bvh4_intersector8_single.h"

#include "geometry/triangle4_intersector8_moeller.h"
#include "geometry/triangle8_intersector8_moeller.h"
#include "geometry/triangle4v_intersector8_pluecker.h"

#define SWITCH_THRESHOLD 5
#define SWITCH_DURING_DOWN_TRAVERSAL 1

namespace embree
{
  namespace isa
  {
    template<int types, bool robust, typename PrimitiveIntersector8>
    void BVH4Intersector8Hybrid<types,robust,PrimitiveIntersector8>::intersect(avxb* valid_i, BVH4* bvh, Ray8& ray)
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
      
      while (1) pop:
      {
        /* pop next node from stack */
        assert(sptr_node > stack_node);
        sptr_node--;
        sptr_near--;
        NodeRef cur = *sptr_node;
        if (unlikely(cur == BVH4::invalidNode)) {
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
            BVH4Intersector8Single<types,robust,PrimitiveIntersector8>::intersect1(bvh, cur, i, pre, ray, ray_org, ray_dir, rdir, ray_tnear, ray_tfar, nearXYZ);
          }
          ray_tfar = min(ray_tfar,ray.tfar);
          continue;
        }
#endif

        while (1)
        {
	  /* process normal nodes */
          if (likely((types & 0x1) && cur.isNode()))
          {
	    const avxb valid_node = ray_tfar > curDist;
	    STAT3(normal.trav_nodes,1,popcnt(valid_node),8);
	    const Node* __restrict__ const node = cur.node();
	    
	    /* pop of next node */
	    assert(sptr_node > stack_node);
	    sptr_node--;
	    sptr_near--;
	    cur = *sptr_node; 
	    curDist = *sptr_near;
	    
#pragma unroll(4)
	    for (unsigned i=0; i<BVH4::N; i++)
	    {
	      const NodeRef child = node->children[i];
	      if (unlikely(child == BVH4::emptyNode)) break;
	      avxf lnearP; const avxb lhit = node->intersect8<robust>(i,org,rdir,org_rdir,ray_tnear,ray_tfar,lnearP);
	      	      
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
#if SWITCH_DURING_DOWN_TRAVERSAL == 1
          // seems to be the best place for testing utilization
          if (unlikely(popcnt(ray_tfar > curDist) <= SWITCH_THRESHOLD))
            {
              *sptr_node++ = cur;
              *sptr_near++ = curDist;
              goto pop;
            }
#endif
	  }
	  
	  /* process motion blur nodes */
          else if (likely((types & 0x10) && cur.isNodeMB()))
	  {
	    const avxb valid_node = ray_tfar > curDist;
	    STAT3(normal.trav_nodes,1,popcnt(valid_node),8);
	    const BVH4::NodeMB* __restrict__ const node = cur.nodeMB();
          
	    /* pop of next node */
	    assert(sptr_node > stack_node);
	    sptr_node--;
	    sptr_near--;
	    cur = *sptr_node; 
	    curDist = *sptr_near;
	    
#pragma unroll(4)
	    for (unsigned i=0; i<BVH4::N; i++)
	    {
	      const NodeRef child = node->child(i);
	      if (unlikely(child == BVH4::emptyNode)) break;
	      avxf lnearP; const avxb lhit = node->intersect(i,org,rdir,org_rdir,ray_tnear,ray_tfar,ray.time,lnearP);
	      	      
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
#if SWITCH_DURING_DOWN_TRAVERSAL == 1
          // seems to be the best place for testing utilization
          if (unlikely(popcnt(ray_tfar > curDist) <= SWITCH_THRESHOLD))
            {
              *sptr_node++ = cur;
              *sptr_near++ = curDist;
              goto pop;
            }
#endif
	  }
	  else 
	    break;
	}
        
        /* return if stack is empty */
        if (unlikely(cur == BVH4::invalidNode)) {
          assert(sptr_node == stack_node);
          break;
        }
        
        /* intersect leaf */
	assert(cur != BVH4::emptyNode);
        const avxb valid_leaf = ray_tfar > curDist;
        STAT3(normal.trav_leaves,1,popcnt(valid_leaf),8);
        size_t items; const Primitive* prim = (Primitive*) cur.leaf(items);
        PrimitiveIntersector8::intersect(valid_leaf,pre,ray,prim,items,bvh->scene);
        ray_tfar = select(valid_leaf,ray.tfar,ray_tfar);
      }
      AVX_ZERO_UPPER();
    }

    
    template<int types, bool robust, typename PrimitiveIntersector8>
    void BVH4Intersector8Hybrid<types,robust,PrimitiveIntersector8>::occluded(avxb* valid_i, BVH4* bvh, Ray8& ray)
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
      
      while (1) pop:
      {
        /* pop next node from stack */
        assert(sptr_node > stack_node);
        sptr_node--;
        sptr_near--;
        NodeRef cur = *sptr_node;
        if (unlikely(cur == BVH4::invalidNode)) {
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
            if (BVH4Intersector8Single<types,robust,PrimitiveIntersector8>::occluded1(bvh,cur,i,pre,ray,ray_org,ray_dir,rdir,ray_tnear,ray_tfar,nearXYZ))
              terminated[i] = -1;
          }
          if (all(terminated)) break;
          ray_tfar = select(terminated,avxf(neg_inf),ray_tfar);
          continue;
        }
#endif
                
        while (1)
        {
	  /* process normal nodes */
          if (likely((types & 0x1) && cur.isNode()))
          {
	    const avxb valid_node = ray_tfar > curDist;
	    STAT3(normal.trav_nodes,1,popcnt(valid_node),8);
	    const Node* __restrict__ const node = cur.node();
	    
	    /* pop of next node */
	    assert(sptr_node > stack_node);
	    sptr_node--;
	    sptr_near--;
	    cur = *sptr_node; 
	    curDist = *sptr_near;
	    
#pragma unroll(4)
	    for (unsigned i=0; i<BVH4::N; i++)
	    {
	      const NodeRef child = node->children[i];
	      if (unlikely(child == BVH4::emptyNode)) break;
	      avxf lnearP; const avxb lhit = node->intersect8<robust>(i,org,rdir,org_rdir,ray_tnear,ray_tfar,lnearP);
	      	      
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
#if SWITCH_DURING_DOWN_TRAVERSAL == 1
          // seems to be the best place for testing utilization
          if (unlikely(popcnt(ray_tfar > curDist) <= SWITCH_THRESHOLD))
            {
              *sptr_node++ = cur;
              *sptr_near++ = curDist;
              goto pop;
            }
#endif
	  }
	  
	  /* process motion blur nodes */
          else if (likely((types & 0x10) && cur.isNodeMB()))
	  {
	    const avxb valid_node = ray_tfar > curDist;
	    STAT3(normal.trav_nodes,1,popcnt(valid_node),8);
	    const BVH4::NodeMB* __restrict__ const node = cur.nodeMB();
          
	    /* pop of next node */
	    assert(sptr_node > stack_node);
	    sptr_node--;
	    sptr_near--;
	    cur = *sptr_node; 
	    curDist = *sptr_near;
	    
#pragma unroll(4)
	    for (unsigned i=0; i<BVH4::N; i++)
	    {
	      const NodeRef child = node->child(i);
	      if (unlikely(child == BVH4::emptyNode)) break;
	      avxf lnearP; const avxb lhit = node->intersect(i,org,rdir,org_rdir,ray_tnear,ray_tfar,ray.time,lnearP);
	      
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
#if SWITCH_DURING_DOWN_TRAVERSAL == 1
          // seems to be the best place for testing utilization
          if (unlikely(popcnt(ray_tfar > curDist) <= SWITCH_THRESHOLD))
            {
              *sptr_node++ = cur;
              *sptr_near++ = curDist;
              goto pop;
            }
#endif
	  }
	  else 
	    break;
	}
        
        /* return if stack is empty */
        if (unlikely(cur == BVH4::invalidNode)) {
          assert(sptr_node == stack_node);
          break;
        }

        
        /* intersect leaf */
	assert(cur != BVH4::emptyNode);
        const avxb valid_leaf = ray_tfar > curDist;
        STAT3(shadow.trav_leaves,1,popcnt(valid_leaf),8);
        size_t items; const Primitive* prim = (Primitive*) cur.leaf(items);
        terminated |= PrimitiveIntersector8::occluded(!terminated,pre,ray,prim,items,bvh->scene);
        if (all(terminated)) break;
        ray_tfar = select(terminated,avxf(neg_inf),ray_tfar);
      }
      store8i(valid & terminated,&ray.geomID,0);
      AVX_ZERO_UPPER();
    }
    
    DEFINE_INTERSECTOR8(BVH4Triangle4Intersector8HybridMoeller, BVH4Intersector8Hybrid<0x1 COMMA false COMMA LeafIterator8_1<Triangle4Intersector8MoellerTrumbore<LeafMode COMMA true> > >);
    DEFINE_INTERSECTOR8(BVH4Triangle4Intersector8HybridMoellerNoFilter, BVH4Intersector8Hybrid<0x1 COMMA false COMMA LeafIterator8_1<Triangle4Intersector8MoellerTrumbore<LeafMode COMMA false> > >);
    DEFINE_INTERSECTOR8(BVH4Triangle8Intersector8HybridMoeller, BVH4Intersector8Hybrid<0x1 COMMA false COMMA LeafIterator8_1<Triangle8Intersector8MoellerTrumbore<LeafMode COMMA true> > >);
    DEFINE_INTERSECTOR8(BVH4Triangle8Intersector8HybridMoellerNoFilter, BVH4Intersector8Hybrid<0x1 COMMA false COMMA LeafIterator8_1<Triangle8Intersector8MoellerTrumbore<LeafMode COMMA false> > >);
    DEFINE_INTERSECTOR8(BVH4Triangle4vIntersector8HybridPluecker, BVH4Intersector8Hybrid<0x1 COMMA true COMMA LeafIterator8_1<Triangle4vIntersector8Pluecker<LeafMode> > >);
  }
}
