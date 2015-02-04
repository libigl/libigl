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

#include "bvh4hair_traversal.h"
#include "bvh4hair_intersector16_single.h"
#include "geometry/bezier1i.h"
#include "geometry/bezier1i_intersector16.h"

#define DBG(x) 

//#define EMBREE_DISABLE_HAIR 0

#if EMBREE_DISABLE_HAIR
// iw, 7/4/14 - Added this workaround to enable external libraries to build
// embree in a way that adds a workardoun for a compiler issue in the
// -mmic branch of icc 13.0 (which breaks on the code that's uncommented below)
#endif

namespace embree
{
  namespace isa
  {
    static unsigned int BVH4HAIR_LEAF_MASK        = BVH4Hair::leaf_mask;        // needed due to compiler efficiency bug
    static unsigned int BVH4HAIR_ALIGNEDNODE_MASK = BVH4Hair::alignednode_mask; // needed due to compiler efficiency bug

  template< bool ENABLE_INTERSECTION_FILTER>
    struct Bezier1iLeafIntersector
    {
      static __forceinline bool intersect(BVH4Hair::NodeRef curNode,
					  const size_t rayIndex, 
					  const mic_f &dir_xyz,
					  const mic_f &org_xyz,
					  const mic_f &min_dist_xyz,
					  mic_f &max_dist_xyz,
					  Ray16& ray16, 
					  const void *__restrict__ const accel,
					  const Scene*__restrict__ const geometry,
					  Precalculations &pre)
      {
	const mic_f pre_vx = broadcast4to16f((float*)&pre.ray_space.vx);
	const mic_f pre_vy = broadcast4to16f((float*)&pre.ray_space.vy);
	const mic_f pre_vz = broadcast4to16f((float*)&pre.ray_space.vz);

	unsigned int items = curNode.items();
	unsigned int index = curNode.offsetIndex();
	const Bezier1i *__restrict__ const tptr = (Bezier1i*)accel + index;
	bool ret = false;
	prefetch<PFHINT_L1>(tptr + 0);
	prefetch<PFHINT_L1>(tptr + 4);

	for (size_t i=0;i<items;i++)
	  ret |= Bezier1iIntersector16<ENABLE_INTERSECTION_FILTER>::intersect(pre_vx,pre_vy,pre_vz,pre.inv_ray_length,ray16,dir_xyz,org_xyz,rayIndex,tptr[i],geometry); 

	max_dist_xyz = ray16.tfar[rayIndex];

	return ret;
      }

      static __forceinline bool occluded(BVH4Hair::NodeRef curNode,
					 const size_t rayIndex, 
					 const mic_f &dir_xyz,
					 const mic_f &org_xyz,
					 const mic_f &min_dist_xyz,
					 const mic_f &max_dist_xyz,
					 const Ray16& ray16, 
					 mic_m &m_terminated,
					 const void *__restrict__ const accel,
					 const Scene*__restrict__ const geometry,
					 Precalculations &pre)
      {
	const mic_f pre_vx = broadcast4to16f((float*)&pre.ray_space.vx);
	const mic_f pre_vy = broadcast4to16f((float*)&pre.ray_space.vy);
	const mic_f pre_vz = broadcast4to16f((float*)&pre.ray_space.vz);

	unsigned int items = curNode.items();
	unsigned int index = curNode.offsetIndex();
	const Bezier1i *__restrict__ const tptr = (Bezier1i*)accel + index;
	prefetch<PFHINT_L1>(tptr + 0);
	prefetch<PFHINT_L1>(tptr + 4);

	for (size_t i=0;i<items;i++)
	  {
	    if (Bezier1iIntersector16<ENABLE_INTERSECTION_FILTER>::occluded(pre_vx,pre_vy,pre_vz,pre.inv_ray_length,ray16,dir_xyz,org_xyz,rayIndex,tptr[i],geometry))
	      return true;
	  }

	return false;
      }

      static __forceinline bool intersect(BVH4Hair::NodeRef curNode,
					  const mic_f &dir_xyz,
					  const mic_f &org_xyz,
					  const mic_f &min_dist_xyz,
					  mic_f &max_dist_xyz,
					  Ray& ray, 
					  const void *__restrict__ const accel,
					  const Scene*__restrict__ const geometry,
					  Precalculations &pre)
      {
	const mic_f pre_vx = broadcast4to16f((float*)&pre.ray_space.vx);
	const mic_f pre_vy = broadcast4to16f((float*)&pre.ray_space.vy);
	const mic_f pre_vz = broadcast4to16f((float*)&pre.ray_space.vz);

	unsigned int items = curNode.items();
	unsigned int index = curNode.offsetIndex(); 
	const Bezier1i *__restrict__ const tptr = (Bezier1i*)accel + index;

	int old_primID = ray.primID;

	bool ret = false;
	prefetch<PFHINT_L1>(tptr + 0);
	prefetch<PFHINT_L1>(tptr + 4);

	for (size_t i=0;i<items;i++)
	  ret |= Bezier1iIntersector16<ENABLE_INTERSECTION_FILTER>::intersect(pre_vx,pre_vy,pre_vz,pre.inv_ray_length,ray,dir_xyz,org_xyz,tptr[i],geometry); 

	max_dist_xyz = ray.tfar;
	return ret;

	return old_primID != ray.primID;
      }

      static __forceinline bool occluded(BVH4Hair::NodeRef curNode,
					 const mic_f &dir_xyz,
					 const mic_f &org_xyz,
					 const mic_f &min_dist_xyz,
					 const mic_f &max_dist_xyz,
					 Ray& ray,
					 const void *__restrict__ const accel,
					 const Scene*__restrict__ const geometry,
					 Precalculations &pre)
      {
	const mic_f pre_vx = broadcast4to16f((float*)&pre.ray_space.vx);
	const mic_f pre_vy = broadcast4to16f((float*)&pre.ray_space.vy);
	const mic_f pre_vz = broadcast4to16f((float*)&pre.ray_space.vz);

	unsigned int items = curNode.items();
	unsigned int index = curNode.offsetIndex(); 
	const Bezier1i *__restrict__ const tptr = (Bezier1i*)accel + index;

	prefetch<PFHINT_L1>(tptr + 0);
	prefetch<PFHINT_L1>(tptr + 4);

	for (size_t i=0;i<items;i++)
	  if (Bezier1iIntersector16<ENABLE_INTERSECTION_FILTER>::occluded(pre_vx,pre_vy,pre_vz,pre.inv_ray_length,ray,dir_xyz,org_xyz,tptr[i],geometry))
	    return true;

	return false;
      }

    };



    template<typename LeafIntersector>    
    void BVH4HairIntersector16<LeafIntersector>::intersect(mic_i* valid_i, BVH4Hair* bvh, Ray16& ray16)
    {
#if EMBREE_DISABLE_HAIR
	THROW_RUNTIME_ERROR("hair explicitly disabled (to work aroudn compiler bug in icc 13.1.0)");
#else 
      /* near and node stack */
      __aligned(64) float   stack_dist[3*BVH4Hair::maxDepth+1];
      __aligned(64) BVH4Hair::NodeRef stack_node[3*BVH4Hair::maxDepth+1];

      const mic_f inv_ray_length = rsqrt(dot(ray16.dir,ray16.dir));
      const mic3f ray16_normalized = ray16.dir * inv_ray_length;
      LinearSpace_mic3f ray16_space = frame(ray16_normalized).transposed();

      /* setup */
      const mic_m m_valid    = *(mic_i*)valid_i != mic_i(0);
      const mic3f rdir16     = rcp_safe(ray16.dir);
      const mic_f inf        = mic_f(pos_inf);
      const mic_f zero       = mic_f::zero();

      store16f(stack_dist,inf);

      const void * __restrict__ accel = (void*)bvh->primitivesPtr();
      const void * __restrict__ nodes = (void*)bvh->nodePtr();

      stack_node[0] = BVH4Hair::invalidNode;
      long rayIndex = -1;
      while((rayIndex = bitscan64(rayIndex,toInt(m_valid))) != BITSCAN_NO_BIT_SET_64)	    
        {
	  Precalculations pre(ray16_space,inv_ray_length,rayIndex);
	  DBG(std::cout << std::endl);
	  DBG(DBG_PRINT(rayIndex));

	  stack_node[1] = bvh->root; 

	  size_t sindex = 2;

	  const mic_f org_xyz      = loadAOS4to16f(rayIndex,ray16.org.x,ray16.org.y,ray16.org.z);
	  const mic_f dir_xyz      = loadAOS4to16f(rayIndex,ray16.dir.x,ray16.dir.y,ray16.dir.z);
	  const mic_f min_dist_xyz = broadcast1to16f(&ray16.tnear[rayIndex]);
	  mic_f       max_dist_xyz = broadcast1to16f(&ray16.tfar[rayIndex]);

	  const mic_f rdir_xyz     = loadAOS4to16f(rayIndex,rdir16.x,rdir16.y,rdir16.z);
	  const mic_f org_rdir_xyz = rdir_xyz * org_xyz;

	  const unsigned int leaf_mask = BVH4HAIR_LEAF_MASK;
	  const unsigned int alignednode_mask = BVH4HAIR_ALIGNEDNODE_MASK;

	  while (1)
	    {

	      BVH4Hair::NodeRef curNode = stack_node[sindex-1];
	      sindex--;

	      traverse_single_intersect(curNode,
					sindex,
					dir_xyz,
					org_xyz,
					rdir_xyz,
					org_rdir_xyz,
					min_dist_xyz,
					max_dist_xyz,
					stack_node,
					stack_dist,
					nodes,
					leaf_mask,
					alignednode_mask);
	      
	      /* return if stack is empty */
	      if (unlikely(curNode == BVH4Hair::invalidNode)) break;

	      STAT3(normal.trav_leaves,1,1,1);

	      /* intersect one ray against bezier curves */

	      //////////////////////////////////////////////////////////////////////////////////////////////////
	      BVH4Hair::NodeRef curNode4i = (unsigned int)curNode;

	      DBG(DBG_PRINT(curNode));

	      const bool hit = LeafIntersector::intersect(curNode4i,
							  rayIndex,
							  dir_xyz,
							  org_xyz,
							  min_dist_xyz,
							  max_dist_xyz,
							  ray16,
							  accel,
							  (Scene*)bvh->geometry,
							  pre);
									   
	      if (hit) 
		{
		  const unsigned int current_dist = *(unsigned int*)&ray16.tfar[rayIndex];
		  compactStack(stack_node,stack_dist,sindex,current_dist,max_dist_xyz);
		}
	      // ------------------------
	    }	  
	}
#endif
    }

    template<typename LeafIntersector>    
    void BVH4HairIntersector16<LeafIntersector>::occluded(mic_i* valid_i, BVH4Hair* bvh, Ray16& ray16)
    {
      /* near and node stack */
      __aligned(64) BVH4Hair::NodeRef stack_node[3*BVH4Hair::maxDepth+1];

      /* setup */
      const mic_f inv_ray_length = rsqrt(dot(ray16.dir,ray16.dir));
      const mic3f ray16_normalized = ray16.dir * inv_ray_length;
      LinearSpace_mic3f ray16_space = frame(ray16_normalized).transposed();

      const mic_m m_valid = *(mic_i*)valid_i != mic_i(0);
      const mic3f rdir16  = rcp_safe(ray16.dir);
      mic_m terminated    = !m_valid;
      const mic_f inf     = mic_f(pos_inf);
      const mic_f zero    = mic_f::zero();

      const void * __restrict__ accel = (void*)bvh->primitivesPtr();
      const void * __restrict__ nodes = (void*)bvh->nodePtr();

      stack_node[0] = BVH4Hair::invalidNode;

      long rayIndex = -1;
      while((rayIndex = bitscan64(rayIndex,toInt(m_valid))) != BITSCAN_NO_BIT_SET_64)	    
        {
	  Precalculations pre(ray16_space,inv_ray_length,rayIndex);

	  stack_node[1] = bvh->root; 
	  size_t sindex = 2;

	  const mic_f org_xyz      = loadAOS4to16f(rayIndex,ray16.org.x,ray16.org.y,ray16.org.z);
	  const mic_f dir_xyz      = loadAOS4to16f(rayIndex,ray16.dir.x,ray16.dir.y,ray16.dir.z);
	  const mic_f min_dist_xyz = broadcast1to16f(&ray16.tnear[rayIndex]);
	  const mic_f max_dist_xyz = broadcast1to16f(&ray16.tfar[rayIndex]);
	  const mic_f rdir_xyz     = loadAOS4to16f(rayIndex,rdir16.x,rdir16.y,rdir16.z);
	  const mic_f org_rdir_xyz = rdir_xyz * org_xyz;

	  const unsigned int leaf_mask    = BVH4HAIR_LEAF_MASK;

	  while (1)
	    {
	      BVH4Hair::NodeRef curNode = stack_node[sindex-1];
	      sindex--;

	      traverse_single_occluded(curNode,
				       sindex,
				       dir_xyz,
				       org_xyz,
				       rdir_xyz,
				       org_rdir_xyz,
				       min_dist_xyz,
				       max_dist_xyz,
				       stack_node,
				       nodes,
				       leaf_mask);


	      /* return if stack is empty */
	      if (unlikely(curNode == BVH4Hair::invalidNode)) break;

	      STAT3(shadow.trav_leaves,1,1,1);

	      /* intersect one ray against bezier curves */

	      //////////////////////////////////////////////////////////////////////////////////////////////////
	      BVH4Hair::NodeRef curNode4i = (unsigned int)curNode;

	      const bool hit = LeafIntersector::occluded(curNode4i,
							 rayIndex,
							 dir_xyz,
							 org_xyz,
							 min_dist_xyz,
							 max_dist_xyz,
							 ray16,
							 terminated,
							 accel,
							 (Scene*)bvh->geometry,
							 pre);

	      if (unlikely(hit)) break;
	      //////////////////////////////////////////////////////////////////////////////////////////////////

	    }


	  if (unlikely(all(toMask(terminated)))) break;
	}


      store16i(m_valid & toMask(terminated),&ray16.geomID,0);
    }
    

    template<typename LeafIntersector>        
    void BVH4HairIntersector1<LeafIntersector>::intersect(BVH4Hair* bvh, Ray& ray)
    {
#if EMBREE_DISABLE_HAIR
	THROW_RUNTIME_ERROR("hair explicitly disabled (to work aroudn compiler bug in icc 13.1.0)");
#else 

      /* near and node stack */
      __aligned(64) float   stack_dist[3*BVH4Hair::maxDepth+1];
      __aligned(64) BVH4Hair::NodeRef stack_node[3*BVH4Hair::maxDepth+1];

      const mic3f ray16_dir            = mic3f(ray.dir.x,ray.dir.y,ray.dir.z);
      const mic_f inv_ray16_length     = rsqrt(dot(ray16_dir,ray16_dir));
      const mic3f ray16_dir_normalized = ray16_dir * inv_ray16_length;
      LinearSpace_mic3f ray16_space    = frame(ray16_dir_normalized).transposed();

      /* setup */
      const mic_f inf        = mic_f(pos_inf);
      const mic_f zero       = mic_f::zero();

      store16f(stack_dist,inf);

      const void * __restrict__ accel = (void*)bvh->primitivesPtr();
      const void * __restrict__ nodes = (void*)bvh->nodePtr();

      stack_node[0] = BVH4Hair::invalidNode;

      Precalculations pre(ray16_space,inv_ray16_length,0);
	  
      stack_node[1] = bvh->root; 
      size_t sindex = 2;

      const mic_f org_xyz      = loadAOS4to16f(ray.org.x,ray.org.y,ray.org.z);
      const mic_f dir_xyz      = loadAOS4to16f(ray.dir.x,ray.dir.y,ray.dir.z);
      const mic_f min_dist_xyz = broadcast1to16f(&ray.tnear);
      mic_f       max_dist_xyz = broadcast1to16f(&ray.tfar);

      const mic_f org_xyz1     = select(0x7777,org_xyz,mic_f::one());
      const mic_f rdir_xyz     = rcp_safe(dir_xyz);
      const mic_f org_rdir_xyz = rdir_xyz * org_xyz;

      const unsigned int leaf_mask = BVH4HAIR_LEAF_MASK;
      const unsigned int alignednode_mask = BVH4HAIR_ALIGNEDNODE_MASK;

      while (1)
	{

	  BVH4Hair::NodeRef curNode = stack_node[sindex-1];
	  sindex--;

	  traverse_single_intersect(curNode,
				    sindex,
				    dir_xyz,
				    org_xyz1,
				    rdir_xyz,
				    org_rdir_xyz,
				    min_dist_xyz,
				    max_dist_xyz,
				    stack_node,
				    stack_dist,
				    nodes,
				    leaf_mask,
				    alignednode_mask);

	  /* return if stack is empty */
	  if (unlikely(curNode == BVH4Hair::invalidNode)) break;

	  STAT3(normal.trav_leaves,1,1,1);
	  STAT3(normal.trav_prims,4,4,4);

	  /* intersect one ray against four bezier curves */

	  //////////////////////////////////////////////////////////////////////////////////////////////////
	  BVH4Hair::NodeRef curNode4i = (unsigned int)curNode;
	  const bool hit = LeafIntersector::intersect(curNode4i,
						      dir_xyz,
						      org_xyz,
						      min_dist_xyz,
						      max_dist_xyz,
						      ray,
						      accel,
						      (Scene*)bvh->geometry,
						      pre);
									   
	  if (hit) 
	    {
	      const unsigned int current_dist = *(unsigned int*)&ray.tfar;
	      compactStack(stack_node,stack_dist,sindex,current_dist,max_dist_xyz);
	    }
	  // ------------------------
	}	         
#endif
    }

    template<typename LeafIntersector>    
    void BVH4HairIntersector1<LeafIntersector>::occluded(BVH4Hair* bvh, Ray& ray)
    {
#if EMBREE_DISABLE_HAIR
	THROW_RUNTIME_ERROR("hair explicitly disabled (to work aroudn compiler bug in icc 13.1.0)");
#else 
      /* near and node stack */
      __aligned(64) float   stack_dist[3*BVH4Hair::maxDepth+1];
      __aligned(64) BVH4Hair::NodeRef stack_node[3*BVH4Hair::maxDepth+1];

      const mic3f ray16_dir            = mic3f(ray.dir.x,ray.dir.y,ray.dir.z);
      const mic_f inv_ray16_length     = rsqrt(dot(ray16_dir,ray16_dir));
      const mic3f ray16_dir_normalized = ray16_dir * inv_ray16_length;
      LinearSpace_mic3f ray16_space    = frame(ray16_dir_normalized).transposed();

      /* setup */
      const mic_f inf        = mic_f(pos_inf);
      const mic_f zero       = mic_f::zero();

      store16f(stack_dist,inf);

      const void * __restrict__ accel = (void*)bvh->primitivesPtr();
      const void * __restrict__ nodes = (void*)bvh->nodePtr();

      stack_node[0] = BVH4Hair::invalidNode;

      Precalculations pre(ray16_space,inv_ray16_length,0);
	  
      stack_node[1] = bvh->root; 
      size_t sindex = 2;

      const mic_f org_xyz      = loadAOS4to16f(ray.org.x,ray.org.y,ray.org.z);
      const mic_f dir_xyz      = loadAOS4to16f(ray.dir.x,ray.dir.y,ray.dir.z);
      const mic_f rdir_xyz     = rcp_safe(dir_xyz);
      const mic_f org_rdir_xyz = rdir_xyz * org_xyz;

      const mic_f min_dist_xyz = broadcast1to16f(&ray.tnear);
      mic_f       max_dist_xyz = broadcast1to16f(&ray.tfar);

      const mic_f org_xyz1     = select(0x7777,org_xyz,mic_f::one());

      const size_t leaf_mask = BVH4HAIR_LEAF_MASK;

      while (1)
	{

	  BVH4Hair::NodeRef curNode = stack_node[sindex-1];
	  sindex--;

	  traverse_single_occluded(curNode,
				   sindex,
				   dir_xyz,
				   org_xyz1,
				   rdir_xyz,
				   org_rdir_xyz,
				   min_dist_xyz,
				   max_dist_xyz,
				   stack_node,
				   nodes,
				   leaf_mask);

	  /* return if stack is empty */
	  if (unlikely(curNode == BVH4Hair::invalidNode)) break;

	  STAT3(normal.trav_leaves,1,1,1);
	  STAT3(normal.trav_prims,4,4,4);

	  /* intersect one ray against bezier curves */

	  //////////////////////////////////////////////////////////////////////////////////////////////////
	  BVH4Hair::NodeRef curNode4i = (unsigned int)curNode;
	  const bool hit = LeafIntersector::occluded(curNode4i,
						     dir_xyz,
						     org_xyz,
						     min_dist_xyz,
						     max_dist_xyz,
						     ray,
						     accel,
						     (Scene*)bvh->geometry,
						     pre);
									   
	  if (unlikely(hit)) 
	    {
	      ray.geomID = 0;
	      return;
	    }
	  // ------------------------
	}	         
#endif
    }
    

    DEFINE_INTERSECTOR1    (BVH4HairIntersector1Bezier1i , BVH4HairIntersector1< Bezier1iLeafIntersector<true> >);
    DEFINE_INTERSECTOR1    (BVH4HairIntersector1Bezier1iNoFilter , BVH4HairIntersector1< Bezier1iLeafIntersector<false> >);
    
    DEFINE_INTERSECTOR16   (BVH4HairIntersector16Bezier1i, BVH4HairIntersector16< Bezier1iLeafIntersector<true> >);
    DEFINE_INTERSECTOR16   (BVH4HairIntersector16Bezier1iNoFilter, BVH4HairIntersector16< Bezier1iLeafIntersector<false> >);

  }
}
