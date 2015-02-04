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

#include "bvh4i.h"
#include "geometry/triangle1.h"
#include "geometry/triangle1_intersector1_moeller.h"
#include "geometry/triangle1mc_intersector1_moeller.h"
#include "geometry/triangle1_intersector16_moeller.h"
#include "geometry/triangle1mc_intersector16_moeller.h"
#include "geometry/virtual_accel_intersector1.h"
#include "geometry/virtual_accel_intersector16.h"
#include "geometry/subdiv_intersector16.h"
#include "geometry/subdiv_intersector1.h"
#include "geometry/filter.h"


namespace embree
{

  // ============================================================================================
  // ============================================================================================
  // ============================================================================================
  static __aligned(64) int zlc4[4] = {0xffffffff,0xffffffff,0xffffffff,0};

  template<bool ENABLE_INTERSECTION_FILTER>
    struct Triangle1LeafIntersector
    {
      // ==================
      // === single ray === 
      // ==================
      static __forceinline bool intersect(BVH4i::NodeRef curNode,
					  const mic_f &dir_xyz,
					  const mic_f &org_xyz,
					  const mic_f &min_dist_xyz,
					  mic_f &max_dist_xyz,
					  Ray& ray, 
					  const void *__restrict__ const accel,
					  const Scene*__restrict__ const geometry)
      {
	const Triangle1* __restrict__ const tptr  = (Triangle1*) curNode.leaf(accel);	      
	const mic_i and_mask = broadcast4to16i(zlc4);
	return Triangle1Intersector1MoellerTrumbore<ENABLE_INTERSECTION_FILTER>::intersect1(dir_xyz,
											    org_xyz,
											    min_dist_xyz,
											    max_dist_xyz,
											    and_mask,
											    ray,
											    geometry,
											    tptr);	
      }

      static __forceinline bool occluded(BVH4i::NodeRef curNode,
					 const mic_f &dir_xyz,
					 const mic_f &org_xyz,
					 const mic_f &min_dist_xyz,
					 const mic_f &max_dist_xyz,
					 Ray& ray,
					 const void *__restrict__ const accel,
					 const Scene*__restrict__ const geometry)
      {
	const Triangle1* __restrict__ const tptr  = (Triangle1*) curNode.leaf(accel);	      
	const mic_i and_mask = broadcast4to16i(zlc4);
	return any(Triangle1Intersector1MoellerTrumbore<ENABLE_INTERSECTION_FILTER>::occluded1(dir_xyz,
											       org_xyz,
											       min_dist_xyz,
											       max_dist_xyz,
											       and_mask,
											       ray,
											       geometry,
											       tptr));	
      }

      // ============================================
      // ==== single ray mode for 16-wide packets ===
      // ============================================
      static __forceinline bool intersect(BVH4i::NodeRef curNode,
					  const size_t rayIndex, 
					  const mic_f &dir_xyz,
					  const mic_f &org_xyz,
					  const mic_f &min_dist_xyz,
					  mic_f &max_dist_xyz,
					  Ray16& ray16, 
					  const void *__restrict__ const accel,
					  const Scene*__restrict__ const geometry)
      {
	const Triangle1* __restrict__ const tptr  = (Triangle1*) curNode.leaf(accel);	      
	const mic_i and_mask = broadcast4to16i(zlc4);
	return Triangle1Intersector16MoellerTrumbore<ENABLE_INTERSECTION_FILTER>::intersect1(curNode,
											     rayIndex,
											     dir_xyz,
											     org_xyz,
											     min_dist_xyz,
											     max_dist_xyz,
											     and_mask,
											     ray16,
											     geometry,
											     tptr);	
      }



      static __forceinline bool occluded(BVH4i::NodeRef curNode,
					 const size_t rayIndex, 
					 const mic_f &dir_xyz,
					 const mic_f &org_xyz,
					 const mic_f &min_dist_xyz,
					 const mic_f &max_dist_xyz,
					 const Ray16& ray16, 
					 mic_m &m_terminated,
					 const void *__restrict__ const accel,
					 const Scene*__restrict__ const geometry)
      {
	const Triangle1* __restrict__ const tptr  = (Triangle1*) curNode.leaf(accel);	      
	const mic_i and_mask = broadcast4to16i(zlc4);
	return Triangle1Intersector16MoellerTrumbore<ENABLE_INTERSECTION_FILTER>::occluded1(curNode,
											    rayIndex,
											    dir_xyz,
											    org_xyz,
											    min_dist_xyz,
											    max_dist_xyz,
											    and_mask,
											    ray16,
											    m_terminated,
											    geometry,
											    tptr);	
      }

      // ========================
      // ==== 16-wide packets ===
      // ========================

      __forceinline static void intersect16(BVH4i::NodeRef curNode,
					    const mic_m m_valid_leaf, 
					    const mic3f &dir,
					    const mic3f &org,
					    Ray16& ray16, 
					    const void *__restrict__ const accel,
					    const Scene     *__restrict__ const geometry)
      {
	unsigned int items; 
	const Triangle1* __restrict__ const tptr  = (Triangle1*) curNode.leaf(accel,items);
	Triangle1Intersector16MoellerTrumbore<ENABLE_INTERSECTION_FILTER>::intersect16(m_valid_leaf,items,dir,org,ray16,geometry,tptr);	
      }

      __forceinline static void occluded16(BVH4i::NodeRef curNode,
					   const mic_m m_valid_leaf, 
					   const mic3f &dir,
					   const mic3f &org,
					   Ray16& ray16, 
					   mic_m &m_terminated,					    
					   const void *__restrict__ const accel,
					   const Scene     *__restrict__ const geometry)
      {
	unsigned int items; 
	const Triangle1* __restrict__ const tptr  = (Triangle1*) curNode.leaf(accel,items);
	Triangle1Intersector16MoellerTrumbore<ENABLE_INTERSECTION_FILTER>::occluded16(m_valid_leaf,items,dir,org,ray16,m_terminated,geometry,tptr);
      }


    };


  // ============================================================================================
  // ============================================================================================
  // ============================================================================================

  template<bool ENABLE_INTERSECTION_FILTER>
    struct Triangle1mcLeafIntersector
    {
      // ==================
      // === single ray === 
      // ==================
      static __forceinline bool intersect(BVH4i::NodeRef curNode,
					  const mic_f &dir_xyz,
					  const mic_f &org_xyz,
					  const mic_f &min_dist_xyz,
					  mic_f &max_dist_xyz,
					  Ray& ray, 
					  const void *__restrict__ const accel,
					  const Scene*__restrict__ const geometry)
      {
	unsigned int items = curNode.items();
	unsigned int index = curNode.offsetIndex();
	const Triangle1mc *__restrict__ const tptr = (Triangle1mc*)accel + index;
	const mic_i and_mask = broadcast4to16i(zlc4);

	return Triangle1mcIntersector1MoellerTrumbore<ENABLE_INTERSECTION_FILTER>::intersect1(dir_xyz,
											      org_xyz,
											      min_dist_xyz,
											      max_dist_xyz,
											      and_mask,
											      ray,
											      geometry,
											      tptr);	
      }

      static __forceinline bool occluded(BVH4i::NodeRef curNode,
					 const mic_f &dir_xyz,
					 const mic_f &org_xyz,
					 const mic_f &min_dist_xyz,
					 const mic_f &max_dist_xyz,
					 Ray& ray,
					 const void *__restrict__ const accel,
					 const Scene*__restrict__ const geometry)
      {
	unsigned int items = curNode.items();
	unsigned int index = curNode.offsetIndex();
	const Triangle1mc *__restrict__ const tptr = (Triangle1mc*)accel + index;
	const mic_i and_mask = broadcast4to16i(zlc4);

	return Triangle1mcIntersector1MoellerTrumbore<ENABLE_INTERSECTION_FILTER>::occluded1(dir_xyz,
											     org_xyz,
											     min_dist_xyz,
											     max_dist_xyz,
											     and_mask,
											     ray,
											     geometry,
											     tptr);	
      }

      // ============================================
      // ==== single ray mode for 16-wide packets ===
      // ============================================
      static __forceinline bool intersect(BVH4i::NodeRef curNode,
					  const size_t rayIndex, 
					  const mic_f &dir_xyz,
					  const mic_f &org_xyz,
					  const mic_f &min_dist_xyz,
					  mic_f &max_dist_xyz,
					  Ray16& ray16, 
					  const void *__restrict__ const accel,
					  const Scene*__restrict__ const geometry)
      {
	unsigned int items = curNode.items();
	unsigned int index = curNode.offsetIndex();
	const Triangle1mc *__restrict__ const tptr = (Triangle1mc*)accel + index;

	const mic_i and_mask = broadcast4to16i(zlc4);
	return Triangle1mcIntersector16MoellerTrumbore<ENABLE_INTERSECTION_FILTER>::intersect1(rayIndex,
											       dir_xyz,
											       org_xyz,
											       min_dist_xyz,
											       max_dist_xyz,
											       and_mask,
											       ray16,
											       geometry,
											       tptr);	
      }


      static __forceinline bool occluded(BVH4i::NodeRef curNode,
					 const size_t rayIndex, 
					 const mic_f &dir_xyz,
					 const mic_f &org_xyz,
					 const mic_f &min_dist_xyz,
					 const mic_f &max_dist_xyz,
					 const Ray16& ray16, 
					 mic_m &m_terminated,
					 const void *__restrict__ const accel,
					 const Scene*__restrict__ const geometry)
      {
	unsigned int items = curNode.items();
	unsigned int index = curNode.offsetIndex();
	const Triangle1mc *__restrict__ const tptr = (Triangle1mc*)accel + index;

	const mic_i and_mask = broadcast4to16i(zlc4);
	return Triangle1mcIntersector16MoellerTrumbore<ENABLE_INTERSECTION_FILTER>::occluded1(rayIndex,
											      dir_xyz,
											      org_xyz,
											      min_dist_xyz,
											      max_dist_xyz,
											      and_mask,
											      ray16,
											      m_terminated,
											      geometry,
											      tptr);	
      }

      // ========================
      // ==== 16-wide packets ===
      // ========================
      __forceinline static void intersect16(BVH4i::NodeRef curNode,
					    const mic_m m_valid_leaf, 
					    const mic3f &dir,
					    const mic3f &org,
					    Ray16& ray16, 
					    const void *__restrict__ const accel,
					    const Scene     *__restrict__ const geometry)
      {
	unsigned int items = curNode.items();
	unsigned int index = curNode.offsetIndex();
	const Triangle1mc *__restrict__ const tptr = (Triangle1mc*)accel + index;
	Triangle1mcIntersector16MoellerTrumbore<ENABLE_INTERSECTION_FILTER>::intersect16(m_valid_leaf,items,dir,org,ray16,geometry,tptr);
      }

      __forceinline static void occluded16(BVH4i::NodeRef curNode,
					   const mic_m m_valid_leaf, 
					   const mic3f &dir,
					   const mic3f &org,
					   Ray16& ray16, 
					   mic_m &m_terminated,					    
					   const void *__restrict__ const accel,
					   const Scene     *__restrict__ const geometry)
      {
	unsigned int items = curNode.items();
	unsigned int index = curNode.offsetIndex();
	const Triangle1mc *__restrict__ const tptr = (Triangle1mc*)accel + index;
	Triangle1mcIntersector16MoellerTrumbore<ENABLE_INTERSECTION_FILTER>::occluded16(m_valid_leaf,items,dir,org,ray16,m_terminated,geometry,tptr);
      }

    };


    // ============================================================================================
    // ============================================================================================
    // ============================================================================================

    template<bool ENABLE_INTERSECTION_FILTER>
    struct VirtualLeafIntersector
    {
      typedef typename VirtualAccelIntersector16::Primitive Primitive;

      // ==================
      // === single ray === 
      // ==================
      static __forceinline bool intersect(BVH4i::NodeRef curNode,
					  const mic_f &dir_xyz,
					  const mic_f &org_xyz,
					  const mic_f &min_dist_xyz,
					  mic_f &max_dist_xyz,
					  Ray& ray, 
					  const void *__restrict__ const accel,
					  const Scene*__restrict__ const geometry)
      {
	unsigned int items = curNode.items();
	unsigned int index = curNode.offsetIndex(); /* array of AccelSetItems */
	Primitive *accel_ptr = (Primitive*)accel + index;
	int old_primID = ray.primID;
	VirtualAccelIntersector1::intersect(ray,accel_ptr,items,geometry);
	return old_primID != ray.primID;
      }

      static __forceinline bool occluded(BVH4i::NodeRef curNode,
					 const mic_f &dir_xyz,
					 const mic_f &org_xyz,
					 const mic_f &min_dist_xyz,
					 const mic_f &max_dist_xyz,
					 Ray& ray,
					 const void *__restrict__ const accel,
					 const Scene*__restrict__ const geometry)
      {
	unsigned int items = curNode.items();
	unsigned int index = curNode.offsetIndex(); /* array of AccelSetItems */
	Primitive *accel_ptr = (Primitive*)accel + index;

	return VirtualAccelIntersector1::occluded((Ray&)ray,accel_ptr,items,(Scene*)geometry);
      }

      // ========================
      // ==== 16-wide packets ===
      // ========================
      __forceinline static void intersect16(BVH4i::NodeRef curNode,
					    const mic_m m_valid_leaf, 
					    const mic3f &dir,
					    const mic3f &org,
					    Ray16& ray16, 
					    const void *__restrict__ const accel,
					    const Scene     *__restrict__ const geometry)
      {
	unsigned int items = curNode.items();
	unsigned int index = curNode.offsetIndex(); /* array of AccelSetItems */
	Primitive *accel_ptr = (Primitive*)accel + index;
        VirtualAccelIntersector16::intersect(m_valid_leaf,ray16,accel_ptr,items,geometry);
      }

      __forceinline static void occluded16(BVH4i::NodeRef curNode,
					   const mic_m m_valid_leaf, 
					   const mic3f &dir,
					   const mic3f &org,
					   Ray16& ray16, 
					   mic_m &m_terminated,					    
					   const void *__restrict__ const accel,
					   const Scene     *__restrict__ const geometry)
      {
	unsigned int items = curNode.items();
	unsigned int index = curNode.offsetIndex(); /* array of AccelSetItems */
	Primitive *accel_ptr = (Primitive *)accel + index;

        m_terminated |= m_valid_leaf & VirtualAccelIntersector16::occluded(m_valid_leaf,ray16,accel_ptr,items,geometry);
      }      
    };





};
