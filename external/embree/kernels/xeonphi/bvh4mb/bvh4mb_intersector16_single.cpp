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

#include "bvh4mb_intersector16_single.h"
#include "bvh4mb_traversal.h"
#include "bvh4mb_leaf_intersector.h"

namespace embree
{
  namespace isa
  {

    static unsigned int BVH4I_LEAF_MASK = BVH4i::leaf_mask; // needed due to compiler efficiency bug

    template<typename LeafIntersector>
    void BVH4mbIntersector16Single<LeafIntersector>::intersect(mic_i* valid_i, BVH4mb* bvh, Ray16& ray16)
    {      
      /* near and node stack */
      __aligned(64) float   stack_dist[3*BVH4i::maxDepth+1];
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

      /* setup */
      const mic_m m_valid    = *(mic_i*)valid_i != mic_i(0);
      const mic3f rdir16     = rcp_safe(ray16.dir);
      const mic_f inf        = mic_f(pos_inf);
      const mic_f zero       = mic_f::zero();
      store16f(stack_dist,inf);

      const Node               * __restrict__ nodes = (Node     *)bvh->nodePtr();
      const BVH4mb::Triangle01 * __restrict__ accel = (BVH4mb::Triangle01 *)bvh->triPtr();

      stack_node[0] = BVH4i::invalidNode;
      long rayIndex = -1;
      while((rayIndex = bitscan64(rayIndex,toInt(m_valid))) != BITSCAN_NO_BIT_SET_64)	    
        {
	  stack_node[1] = bvh->root;
	  size_t sindex = 2;

	  const mic_f org_xyz      = loadAOS4to16f(rayIndex,ray16.org.x,ray16.org.y,ray16.org.z);
	  const mic_f dir_xyz      = loadAOS4to16f(rayIndex,ray16.dir.x,ray16.dir.y,ray16.dir.z);
	  const mic_f rdir_xyz     = loadAOS4to16f(rayIndex,rdir16.x,rdir16.y,rdir16.z);
	  const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
	  const mic_f min_dist_xyz = broadcast1to16f(&ray16.tnear[rayIndex]);
	  mic_f       max_dist_xyz = broadcast1to16f(&ray16.tfar[rayIndex]);
	  const mic_f time         = broadcast1to16f(&ray16.time[rayIndex]);
	      
	  
	  const unsigned int leaf_mask = BVH4I_LEAF_MASK;

	  while (1)
	    {
	      NodeRef curNode = stack_node[sindex-1];
	      sindex--;

	      traverse_single_intersect(curNode,
					sindex,
					rdir_xyz,
					org_rdir_xyz,
					min_dist_xyz,
					max_dist_xyz,
					time,
					stack_node,
					stack_dist,
					nodes,
					leaf_mask);

	    
	      /* return if stack is empty */
	      if (unlikely(curNode == BVH4i::invalidNode)) break;


	      /* intersect one ray against four triangles */

	      //////////////////////////////////////////////////////////////////////////////////////////////////

	      const bool hit = LeafIntersector::intersect(curNode,
							  rayIndex,
							  dir_xyz,
							  org_xyz,
							  min_dist_xyz,
							  max_dist_xyz,
							  ray16,
							  accel,
							  (Scene*)bvh->geometry);
									   
	      if (hit)
		compactStack(stack_node,stack_dist,sindex,max_dist_xyz);

	    }	  
	}
    }
    
    template<typename LeafIntersector>
    void BVH4mbIntersector16Single<LeafIntersector>::occluded(mic_i* valid_i, BVH4mb* bvh, Ray16& ray16)
    {
      /* near and node stack */
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

      /* setup */
      const mic_m m_valid     = *(mic_i*)valid_i != mic_i(0);
      const mic3f rdir16      = rcp_safe(ray16.dir);
      mic_m m_terminated      = !m_valid;
      const mic_f inf         = mic_f(pos_inf);
      const mic_f zero        = mic_f::zero();

      const Node               * __restrict__ nodes = (Node     *)bvh->nodePtr();
      const BVH4mb::Triangle01 * __restrict__ accel = (BVH4mb::Triangle01 *)bvh->triPtr();

      stack_node[0] = BVH4i::invalidNode;

      long rayIndex = -1;
      while((rayIndex = bitscan64(rayIndex,toInt(m_valid))) != BITSCAN_NO_BIT_SET_64)	    
        {
	  stack_node[1] = bvh->root;
	  size_t sindex = 2;

	  const mic_f org_xyz      = loadAOS4to16f(rayIndex,ray16.org.x,ray16.org.y,ray16.org.z);
	  const mic_f dir_xyz      = loadAOS4to16f(rayIndex,ray16.dir.x,ray16.dir.y,ray16.dir.z);
	  const mic_f rdir_xyz     = loadAOS4to16f(rayIndex,rdir16.x,rdir16.y,rdir16.z);
	  const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
	  const mic_f min_dist_xyz = broadcast1to16f(&ray16.tnear[rayIndex]);
	  const mic_f max_dist_xyz = broadcast1to16f(&ray16.tfar[rayIndex]);
	  const mic_f time         = broadcast1to16f(&ray16.time[rayIndex]);

	  const unsigned int leaf_mask = BVH4I_LEAF_MASK;
	  while (1)
	    {
	      NodeRef curNode = stack_node[sindex-1];
	      sindex--;

	      traverse_single_occluded(curNode,
				       sindex,
				       rdir_xyz,
				       org_rdir_xyz,
				       min_dist_xyz,
				       max_dist_xyz,
				       time,
				       stack_node,
				       nodes,
				       leaf_mask);

	    
	      /* return if stack is empty */
	      if (unlikely(curNode == BVH4i::invalidNode)) break;


	      /* intersect one ray against four triangles */

	      //////////////////////////////////////////////////////////////////////////////////////////////////

	      const bool hit = LeafIntersector::occluded(curNode,
							 rayIndex,
							 dir_xyz,
							 org_xyz,
							 min_dist_xyz,
							 max_dist_xyz,
							 ray16,
							 m_terminated,
							 accel,
							 (Scene*)bvh->geometry);

	      if (unlikely(hit)) break;
	    }

	  if (unlikely(all(toMask(m_terminated)))) break;
	}


      store16i(m_valid & m_terminated,&ray16.geomID,0);
    }
    
    DEFINE_INTERSECTOR16    (BVH4mbTriangle1Intersector16SingleMoeller, BVH4mbIntersector16Single<Triangle1mbLeafIntersector>);

  }
}
