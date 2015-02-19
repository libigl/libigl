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

#include "bvh4i_intersector16_chunk.h"
#include "bvh4i_leaf_intersector.h"

namespace embree
{
  namespace isa
  {
    static unsigned int BVH4I_LEAF_MASK = BVH4i::leaf_mask; // needed due to compiler efficiency bug

    template<typename LeafIntersector, bool ENABLE_COMPRESSED_BVH4I_NODES>
    void BVH4iIntersector16Chunk<LeafIntersector,ENABLE_COMPRESSED_BVH4I_NODES>::intersect(mic_i* valid_i, BVH4i* bvh, Ray16& ray)
    {
      /* near and node stack */
      __aligned(64) mic_f   stack_dist[3*BVH4i::maxDepth+1];
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

      /* load ray */
      const mic_m valid0   = *(mic_i*)valid_i != mic_i(0);
      const mic3f rdir     = rcp_safe(ray.dir);
 
      const mic3f org_rdir = ray.org * rdir;
      mic_f ray_tnear      = select(valid0,ray.tnear,pos_inf);
      mic_f ray_tfar       = select(valid0,ray.tfar ,neg_inf);
      const mic_f inf      = mic_f(pos_inf);
      
      /* allocate stack and push root node */
      stack_node[0] = BVH4i::invalidNode;
      stack_dist[0] = inf;
      stack_node[1] = bvh->root;
      stack_dist[1] = ray_tnear; 
      NodeRef* __restrict__ sptr_node = stack_node + 2;
      mic_f*   __restrict__ sptr_dist = stack_dist + 2;
      
      const Node      * __restrict__ nodes = (Node     *)bvh->nodePtr();
      const Triangle1 * __restrict__ accel = (Triangle1*)bvh->triPtr();

      const mic3f org = ray.org;
      const mic3f dir = ray.dir;

      while (1)
      {
        /* pop next node from stack */
        NodeRef curNode = *(sptr_node-1);
        mic_f curDist   = *(sptr_dist-1);
        sptr_node--;
        sptr_dist--;
	const mic_m m_stackDist = ray_tfar > curDist;

	/* stack emppty ? */
        if (unlikely(curNode == BVH4i::invalidNode))  break;
        
        /* cull node if behind closest hit point */
        if (unlikely(none(m_stackDist))) {continue;}
	        
	const unsigned int leaf_mask = BVH4I_LEAF_MASK;

	traverse_chunk_intersect<ENABLE_COMPRESSED_BVH4I_NODES>(curNode,
								curDist,
								rdir,
								org_rdir,
								ray_tnear,
								ray_tfar,
								sptr_node,
								sptr_dist,
								nodes,
								leaf_mask);            		    	
        
        /* return if stack is empty */
        if (unlikely(curNode == BVH4i::invalidNode)) break;
        
        /* intersect leaf */
        const mic_m m_valid_leaf = ray_tfar > curDist;
        STAT3(normal.trav_leaves,1,popcnt(m_valid_leaf),16);
 
	LeafIntersector::intersect16(curNode,m_valid_leaf,dir,org,ray,accel,(Scene*)bvh->geometry);

        ray_tfar = select(m_valid_leaf,ray.tfar,ray_tfar);
      }
    }

    template<typename LeafIntersector,bool ENABLE_COMPRESSED_BVH4I_NODES>
    void BVH4iIntersector16Chunk<LeafIntersector,ENABLE_COMPRESSED_BVH4I_NODES>::occluded(mic_i* valid_i, BVH4i* bvh, Ray16& ray)
    {
      /* allocate stack */
      __aligned(64) mic_f    stack_dist[3*BVH4i::maxDepth+1];
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

      /* load ray */
      const mic_m valid = *(mic_i*)valid_i != mic_i(0);
      mic_m m_terminated = !valid;
      const mic3f rdir = rcp_safe(ray.dir);
      const mic3f org_rdir = ray.org * rdir;
      mic_f ray_tnear = select(valid,ray.tnear,pos_inf);
      mic_f ray_tfar  = select(valid,ray.tfar ,neg_inf);
      const mic_f inf = mic_f(pos_inf);
      
      /* push root node */
      stack_node[0] = BVH4i::invalidNode;
      stack_dist[0] = inf;
      stack_node[1] = bvh->root;
      stack_dist[1] = ray_tnear; 
      NodeRef* __restrict__ sptr_node = stack_node + 2;
      mic_f*   __restrict__ sptr_dist = stack_dist + 2;
      
      const Node      * __restrict__ nodes = (Node     *)bvh->nodePtr();
      const Triangle1 * __restrict__ accel = (Triangle1*)bvh->triPtr();

      const mic3f org = ray.org;
      const mic3f dir = ray.dir;

      while (1)
      {
	const mic_m m_active = !m_terminated;

        /* pop next node from stack */
        NodeRef curNode = *(sptr_node-1);
        mic_f curDist   = *(sptr_dist-1);
        sptr_node--;
        sptr_dist--;
	const mic_m m_stackDist = gt(m_active,ray_tfar,curDist);

	/* stack emppty ? */
        if (unlikely(curNode == BVH4i::invalidNode))  break;
        
        /* cull node if behind closest hit point */

        if (unlikely(none(m_stackDist))) { continue; }

	const unsigned int leaf_mask = BVH4I_LEAF_MASK;

	traverse_chunk_occluded<ENABLE_COMPRESSED_BVH4I_NODES>(curNode,
							       curDist,
							       rdir,
							       org_rdir,
							       ray_tnear,
							       ray_tfar,
							       m_active,
							       sptr_node,
							       sptr_dist,
							       nodes,
							       leaf_mask);            		    	
        
        /* return if stack is empty */
        if (unlikely(curNode == BVH4i::invalidNode)) break;
        
        /* intersect leaf */
        mic_m m_valid_leaf = gt(m_active,ray_tfar,curDist);
        STAT3(shadow.trav_leaves,1,popcnt(m_valid_leaf),16);

	LeafIntersector::occluded16(curNode,m_valid_leaf,dir,org,ray,m_terminated,accel,(Scene*)bvh->geometry);

        if (unlikely(all(m_terminated))) break;
        ray_tfar = select(m_terminated,neg_inf,ray_tfar);
      }
      store16i(valid & m_terminated,&ray.geomID,0);
    }

    typedef BVH4iIntersector16Chunk< Triangle1LeafIntersector  < true  >, false > Triangle1Intersector16ChunkMoellerFilter;
    typedef BVH4iIntersector16Chunk< Triangle1LeafIntersector  < false >, false > Triangle1Intersector16ChunkMoellerNoFilter;
    typedef BVH4iIntersector16Chunk< Triangle1mcLeafIntersector< true  >, true  > Triangle1mcIntersector16ChunkMoellerFilter;
    typedef BVH4iIntersector16Chunk< Triangle1mcLeafIntersector< false >, true  > Triangle1mcIntersector16ChunkMoellerNoFilter;
    typedef BVH4iIntersector16Chunk< VirtualLeafIntersector    < true  >, false > VirtualIntersector16ChunkMoellerFilter;
    typedef BVH4iIntersector16Chunk< VirtualLeafIntersector    < false >, false > VirtualIntersector16ChunkMoellerNoFilter;

    DEFINE_INTERSECTOR16    (BVH4iTriangle1Intersector16ChunkMoeller          , Triangle1Intersector16ChunkMoellerFilter);
    DEFINE_INTERSECTOR16    (BVH4iTriangle1Intersector16ChunkMoellerNoFilter  , Triangle1Intersector16ChunkMoellerNoFilter);
    DEFINE_INTERSECTOR16    (BVH4iTriangle1mcIntersector16ChunkMoeller        , Triangle1mcIntersector16ChunkMoellerFilter);
    DEFINE_INTERSECTOR16    (BVH4iTriangle1mcIntersector16ChunkMoellerNoFilter, Triangle1mcIntersector16ChunkMoellerNoFilter);
    DEFINE_INTERSECTOR16   (BVH4iVirtualGeometryIntersector16        , VirtualIntersector16ChunkMoellerFilter);
    DEFINE_INTERSECTOR16   (BVH4iVirtualGeometryIntersector16NoFilter, VirtualIntersector16ChunkMoellerNoFilter);


  }
}
