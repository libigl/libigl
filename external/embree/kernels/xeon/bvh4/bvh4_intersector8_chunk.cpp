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

#include "bvh4_intersector8_chunk.h"

#include "geometry/bezier1v_intersector8.h"
#include "geometry/bezier1i_intersector8.h"
#include "geometry/triangle1_intersector8_moeller.h"
#include "geometry/triangle4_intersector8_moeller.h"
#include "geometry/triangle8_intersector8_moeller.h"
#include "geometry/triangle1v_intersector8_pluecker.h"
#include "geometry/triangle4v_intersector8_pluecker.h"
#include "geometry/triangle4i_intersector8.h"
#include "geometry/virtual_accel_intersector8.h"
#include "geometry/triangle1v_intersector8_moeller_mb.h"
#include "geometry/triangle4v_intersector8_moeller_mb.h"

namespace embree
{
  namespace isa
  {
    template<int types, bool robust, typename PrimitiveIntersector8>
    void BVH4Intersector8Chunk<types, robust, PrimitiveIntersector8>::intersect(avxb* valid_i, BVH4* bvh, Ray8& ray)
    {
      /* load ray */
      const avxb valid0 = *valid_i;
      const avx3f rdir = rcp_safe(ray.dir);
      const avx3f org(ray.org), org_rdir = org * rdir;
      avxf ray_tnear = select(valid0,ray.tnear,pos_inf);
      avxf ray_tfar  = select(valid0,ray.tfar ,neg_inf);
      const avxf inf = avxf(pos_inf);
      Precalculations pre(valid0,ray);

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
        NodeRef cur = *sptr_node;
        if (unlikely(cur == BVH4::invalidNode)) {
          assert(sptr_node == stack_node);
          break;
        }
        
        /* cull node if behind closest hit point */
        avxf curDist = *sptr_near;
        if (unlikely(none(ray_tfar > curDist))) 
          continue;
        
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
    void BVH4Intersector8Chunk<types, robust, PrimitiveIntersector8>::occluded(avxb* valid_i, BVH4* bvh, Ray8& ray)
    {
      /* load ray */
      const avxb valid = *valid_i;
      avxb terminated = !valid;
      const avx3f rdir = rcp_safe(ray.dir);
      const avx3f org(ray.org), org_rdir = org * rdir;
      avxf ray_tnear = select(valid,ray.tnear,pos_inf);
      avxf ray_tfar  = select(valid,ray.tfar ,neg_inf);
      const avxf inf = avxf(pos_inf);
      Precalculations pre(valid,ray);

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
        NodeRef cur = *sptr_node;
        if (unlikely(cur == BVH4::invalidNode)) {
          assert(sptr_node == stack_node);
          break;
        }
        
        /* cull node if behind closest hit point */
        avxf curDist = *sptr_near;
        if (unlikely(none(ray_tfar > curDist))) 
          continue;
        
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
        terminated |= valid_leaf & PrimitiveIntersector8::occluded(valid_leaf,pre,ray,prim,items,bvh->scene);
        if (all(terminated)) break;
        ray_tfar = select(terminated,neg_inf,ray_tfar);
      }
      store8i(valid & terminated,&ray.geomID,0);
      AVX_ZERO_UPPER();
    }

    DEFINE_INTERSECTOR8(BVH4Bezier1vIntersector8Chunk, BVH4Intersector8Chunk<0x1 COMMA false COMMA LeafIterator8<Bezier1vIntersector8<LeafMode> > >);
    DEFINE_INTERSECTOR8(BVH4Bezier1iIntersector8Chunk, BVH4Intersector8Chunk<0x1 COMMA false COMMA LeafIterator8<Bezier1iIntersector8<LeafMode> > >);
    DEFINE_INTERSECTOR8(BVH4Triangle1Intersector8ChunkMoeller, BVH4Intersector8Chunk<0x1 COMMA false COMMA LeafIterator8<Triangle1Intersector8MoellerTrumbore<LeafMode> > >);
    DEFINE_INTERSECTOR8(BVH4Triangle4Intersector8ChunkMoeller, BVH4Intersector8Chunk<0x1 COMMA false COMMA LeafIterator8<Triangle4Intersector8MoellerTrumbore<LeafMode COMMA true> > >);
    DEFINE_INTERSECTOR8(BVH4Triangle4Intersector8ChunkMoellerNoFilter, BVH4Intersector8Chunk<0x1 COMMA false COMMA LeafIterator8<Triangle4Intersector8MoellerTrumbore<LeafMode COMMA false> > >);
    DEFINE_INTERSECTOR8(BVH4Triangle8Intersector8ChunkMoeller, BVH4Intersector8Chunk<0x1 COMMA false COMMA LeafIterator8<Triangle8Intersector8MoellerTrumbore<LeafMode COMMA true> > >);
    DEFINE_INTERSECTOR8(BVH4Triangle8Intersector8ChunkMoellerNoFilter, BVH4Intersector8Chunk<0x1 COMMA false COMMA LeafIterator8<Triangle8Intersector8MoellerTrumbore<LeafMode COMMA false> > >);
    DEFINE_INTERSECTOR8(BVH4Triangle1vIntersector8ChunkPluecker, BVH4Intersector8Chunk<0x1 COMMA true COMMA LeafIterator8<Triangle1vIntersector8Pluecker<LeafMode> > >);
    DEFINE_INTERSECTOR8(BVH4Triangle4vIntersector8ChunkPluecker, BVH4Intersector8Chunk<0x1 COMMA true COMMA LeafIterator8<Triangle4vIntersector8Pluecker<LeafMode> > >);
    DEFINE_INTERSECTOR8(BVH4Triangle4iIntersector8ChunkPluecker, BVH4Intersector8Chunk<0x1 COMMA true COMMA LeafIterator8<Triangle4iIntersector8Pluecker<LeafMode> > >);
    DEFINE_INTERSECTOR8(BVH4VirtualIntersector8Chunk, BVH4Intersector8Chunk<0x1 COMMA false COMMA LeafIterator8<VirtualAccelIntersector8> >);

    DEFINE_INTERSECTOR8(BVH4Triangle1vMBIntersector8ChunkMoeller, BVH4Intersector8Chunk<0x10 COMMA false COMMA LeafIterator8<Triangle1vIntersector8MoellerTrumboreMB<LeafMode> > >);
    DEFINE_INTERSECTOR8(BVH4Triangle4vMBIntersector8ChunkMoeller, BVH4Intersector8Chunk<0x10 COMMA false COMMA LeafIterator8<Triangle4vMBIntersector8MoellerTrumbore<LeafMode COMMA true> > >);
  }
}
