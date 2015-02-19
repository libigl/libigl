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

#include "triangle1.h"
#include "common/ray.h"
#include "geometry/filter.h"

namespace embree
{
  /*! Intersector for individual precomputed triangles with one
   *  ray. This intersector implements a modified version of the
   *  Moeller Trumbore intersector from the paper "Fast, Minimum
   *  Storage Ray-Triangle Intersection". In contrast to the paper we
   *  precalculate some factors and factor the calculations
   *  differently to allow precalculating the cross product e1 x
   *  e2. */

  template< bool ENABLE_INTERSECTION_FILTER>
    struct Triangle1Intersector1MoellerTrumbore
    {

      __forceinline static bool intersect1(const mic_f &dir_xyz,
					   const mic_f &org_xyz,
					   const mic_f &min_dist_xyz,
					   mic_f &max_dist_xyz,
					   const mic_i &and_mask,
					   Ray& ray, 
					   const Scene     *__restrict__ const geometry,
					   const Triangle1 * __restrict__ const tptr)
      {
	const mic_f zero = mic_f::zero();

	prefetch<PFHINT_L1>(tptr + 3);
	prefetch<PFHINT_L1>(tptr + 2);
	prefetch<PFHINT_L1>(tptr + 1);
	prefetch<PFHINT_L1>(tptr + 0); 
	      
	const mic_f v0 = gather_4f_zlc(and_mask,
				       (float*)&tptr[0].v0,
				       (float*)&tptr[1].v0,
				       (float*)&tptr[2].v0,
				       (float*)&tptr[3].v0);
	      
	const mic_f v1 = gather_4f_zlc(and_mask,
				       (float*)&tptr[0].v1,
				       (float*)&tptr[1].v1,
				       (float*)&tptr[2].v1,
				       (float*)&tptr[3].v1);
	      
	const mic_f v2 = gather_4f_zlc(and_mask,
				       (float*)&tptr[0].v2,
				       (float*)&tptr[1].v2,
				       (float*)&tptr[2].v2,
				       (float*)&tptr[3].v2);

	const mic_f e1 = v1 - v0;
	const mic_f e2 = v0 - v2;	     
	const mic_f normal = lcross_zxy(e1,e2);
	const mic_f org = v0 - org_xyz;
	const mic_f odzxy = msubr231(org * swizzle(dir_xyz,_MM_SWIZ_REG_DACB), dir_xyz, swizzle(org,_MM_SWIZ_REG_DACB));
	const mic_f den = ldot3_zxy(dir_xyz,normal);	      
	const mic_f rcp_den = rcp(den);
	const mic_f uu = ldot3_zxy(e2,odzxy); 
	const mic_f vv = ldot3_zxy(e1,odzxy); 
	const mic_f u = uu * rcp_den;
	const mic_f v = vv * rcp_den;

#if defined(RTCORE_BACKFACE_CULLING)
	const mic_m m_init = (mic_m)0x1111 & (den > zero);
#else
	const mic_m m_init = 0x1111;
#endif

	const mic_m valid_u = ge(m_init,u,zero);
	const mic_m valid_v = ge(valid_u,v,zero);
	const mic_m m_aperture = le(valid_v,u+v,mic_f::one()); 

	const mic_f nom = ldot3_zxy(org,normal);

	if (unlikely(none(m_aperture))) return false;
	const mic_f t = rcp_den*nom;

	mic_m m_final  = lt(lt(m_aperture,min_dist_xyz,t),t,max_dist_xyz);

#if defined(RTCORE_RAY_MASK)
	const mic_i rayMask(ray.mask);
	const mic_i triMask = getTriMasks(tptr);
	const mic_m m_ray_mask = (rayMask & triMask) != mic_i::zero();
	m_final &= m_ray_mask;	      
#endif


	//////////////////////////////////////////////////////////////////////////////////////////////////

	/* did the ray hot one of the four triangles? */
	if (unlikely(any(m_final)))
	  {
	    /* intersection filter test */
	    if (ENABLE_INTERSECTION_FILTER) 
	      {
		mic_f org_max_dist_xyz = max_dist_xyz;

		/* did the ray hit one of the four triangles? */
		while (any(m_final)) 
		  {
		    max_dist_xyz  = select(m_final,t,org_max_dist_xyz);
		    const mic_f min_dist = vreduce_min(max_dist_xyz);
		    const mic_m m_dist = eq(min_dist,max_dist_xyz);
		    const size_t vecIndex = bitscan(toInt(m_dist));
		    const size_t triIndex = vecIndex >> 2;
		    const Triangle1  *__restrict__ tri_ptr = tptr + triIndex;
		    const mic_m m_tri = m_dist^(m_dist & (mic_m)((unsigned int)m_dist - 1));
		    const mic_f gnormalx = mic_f(tri_ptr->Ng.x);
		    const mic_f gnormaly = mic_f(tri_ptr->Ng.y);
		    const mic_f gnormalz = mic_f(tri_ptr->Ng.z);
		    const int geomID = tri_ptr->geomID();
		    const int primID = tri_ptr->primID();
                
		    const Geometry* geom = geometry->get(geomID);
		    if (likely(!geom->hasIntersectionFilter1())) 
		      {
			ray.update(m_tri,min_dist,u,v,gnormalx,gnormaly,gnormalz,geomID,primID);
			break;
		      }
                
		    if (runIntersectionFilter1(geom,ray,u,v,min_dist,gnormalx,gnormaly,gnormalz,m_tri,geomID,primID)) {
		      break;
		    }
		    m_final ^= m_tri;
		  }
		max_dist_xyz = ray.tfar;
	      }
	    else
	      {
		max_dist_xyz  = select(m_final,t,max_dist_xyz);
		const mic_f min_dist = vreduce_min(max_dist_xyz);
		const mic_m m_dist = eq(min_dist,max_dist_xyz);

		prefetch<PFHINT_L1EX>((mic_f*)&ray + 0);
		prefetch<PFHINT_L1EX>((mic_f*)&ray + 1);

		const size_t vecIndex = bitscan(toInt(m_dist));
		const size_t triIndex = vecIndex >> 2;

		const Triangle1  *__restrict__ tri_ptr = tptr + triIndex;

		const mic_m m_tri = m_dist^(m_dist & (mic_m)((unsigned int)m_dist - 1));

		const mic_f gnormalx = mic_f(tri_ptr->Ng.x);
		const mic_f gnormaly = mic_f(tri_ptr->Ng.y);
		const mic_f gnormalz = mic_f(tri_ptr->Ng.z);

		max_dist_xyz = min_dist;

		ray.update(m_tri,min_dist,u,v,gnormalx,gnormaly,gnormalz,tri_ptr->geomID(),tri_ptr->primID());
	      }
	  return true;

	}
      return false;
    }


    __forceinline static mic_m occluded1(const mic_f &dir_xyz,
					 const mic_f &org_xyz,
					 const mic_f &min_dist_xyz,
					 const mic_f &max_dist_xyz,
					 const mic_i &and_mask,
					 const Ray& ray, 
					 const Scene     *__restrict__ const geometry,
					 const Triangle1 * __restrict__ const tptr)
    {
      const mic_f zero = mic_f::zero();

      prefetch<PFHINT_L1>(tptr + 3);
      prefetch<PFHINT_L1>(tptr + 2);
      prefetch<PFHINT_L1>(tptr + 1);
      prefetch<PFHINT_L1>(tptr + 0); 

	      
      const mic_f v0 = gather_4f_zlc(and_mask,
				     (float*)&tptr[0].v0,
				     (float*)&tptr[1].v0,
				     (float*)&tptr[2].v0,
				     (float*)&tptr[3].v0);
	      
      const mic_f v1 = gather_4f_zlc(and_mask,
				     (float*)&tptr[0].v1,
				     (float*)&tptr[1].v1,
				     (float*)&tptr[2].v1,
				     (float*)&tptr[3].v1);
	      
      const mic_f v2 = gather_4f_zlc(and_mask,
				     (float*)&tptr[0].v2,
				     (float*)&tptr[1].v2,
				     (float*)&tptr[2].v2,
				     (float*)&tptr[3].v2);

      const mic_f e1 = v1 - v0;
      const mic_f e2 = v0 - v2;	     
      const mic_f normal = lcross_zxy(e1,e2);
      const mic_f org = v0 - org_xyz;
      const mic_f odzxy = msubr231(org * swizzle(dir_xyz,_MM_SWIZ_REG_DACB), dir_xyz, swizzle(org,_MM_SWIZ_REG_DACB));
      const mic_f den = ldot3_zxy(dir_xyz,normal);	      
      const mic_f rcp_den = rcp(den);
      const mic_f uu = ldot3_zxy(e2,odzxy); 
      const mic_f vv = ldot3_zxy(e1,odzxy); 
      const mic_f u = uu * rcp_den;
      const mic_f v = vv * rcp_den;

#if defined(RTCORE_BACKFACE_CULLING)
      const mic_m m_init = (mic_m)0x1111 & (den > zero);
#else
      const mic_m m_init = 0x1111;
#endif

      const mic_m valid_u = ge((mic_m)m_init,u,zero);
      const mic_m valid_v = ge(valid_u,v,zero);
      const mic_m m_aperture = le(valid_v,u+v,mic_f::one()); 

      const mic_f nom = ldot3_zxy(org,normal);
      const mic_f t = rcp_den*nom;
      if (unlikely(none(m_aperture))) return m_aperture;

      mic_m m_final  = lt(lt(m_aperture,min_dist_xyz,t),t,max_dist_xyz);

#if defined(RTCORE_RAY_MASK)
      const mic_i rayMask(ray.mask);
      const mic_i triMask = getTriMasks(tptr); 
      const mic_m m_ray_mask = (rayMask & triMask) != mic_i::zero();
      m_final &= m_ray_mask;	      
#endif

      if (ENABLE_INTERSECTION_FILTER)
	{
	  /* did the ray hit one of the four triangles? */
	  while (any(m_final)) 
	    {
	      const mic_f temp_t  = select(m_final,t,max_dist_xyz);
	      const mic_f min_dist = vreduce_min(temp_t);
	      const mic_m m_dist = eq(min_dist,temp_t);
	      const size_t vecIndex = bitscan(toInt(m_dist));
	      const size_t triIndex = vecIndex >> 2;
	      const Triangle1  *__restrict__ tri_ptr = tptr + triIndex;
	      const mic_m m_tri = m_dist^(m_dist & (mic_m)((unsigned int)m_dist - 1));
	      const mic_f gnormalx = mic_f(tri_ptr->Ng.x);
	      const mic_f gnormaly = mic_f(tri_ptr->Ng.y);
	      const mic_f gnormalz = mic_f(tri_ptr->Ng.z);
	      const int geomID = tri_ptr->geomID();
	      const int primID = tri_ptr->primID();                
	      const Geometry* geom = geometry->get(geomID);

	      if (likely(!geom->hasOcclusionFilter1())) break;
                
	      if (runOcclusionFilter1(geom,(Ray&)ray,u,v,min_dist,gnormalx,gnormaly,gnormalz,m_tri,geomID,primID)) 
		break;

	      m_final ^= m_tri; /* clear bit */
	    }
	}

      return m_final;
    }

  };
}
