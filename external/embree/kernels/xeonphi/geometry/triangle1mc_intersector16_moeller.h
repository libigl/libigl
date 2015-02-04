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
#include "common/ray16.h"
#include "geometry/filter.h"

namespace embree
{
  /*! Intersector for individual precomputed triangles with 16
   *  rays. This intersector implements a modified version of the
   *  Moeller Trumbore intersector from the paper "Fast, Minimum
   *  Storage Ray-Triangle Intersection". In contrast to the paper we
   *  precalculate some factors and factor the calculations
   *  differently to allow precalculating the cross product e1 x
   *  e2. */

  template< bool ENABLE_INTERSECTION_FILTER>
    struct Triangle1mcIntersector16MoellerTrumbore
    {

      __forceinline static bool intersect1(const size_t rayIndex, 
					   const mic_f &dir_xyz,
					   const mic_f &org_xyz,
					   const mic_f &min_dist_xyz,
					   mic_f &max_dist_xyz,
					   const mic_i &and_mask,
					   Ray16& ray16, 
					   const Scene     *__restrict__ const scene,
					   const Triangle1mc * __restrict__ const tptr)
      {
	const mic_f zero = mic_f::zero();
	prefetch<PFHINT_L1>(tptr + 0); 
	prefetch<PFHINT_L1>(tptr + 2);


	const mic_f v0 = gather_4f_zlc(and_mask,
				       tptr[0].v0,
				       tptr[1].v0,
				       tptr[2].v0,
				       tptr[3].v0);
	      
	const mic_f v1 = gather_4f_zlc(and_mask,
				       tptr[0].v1,
				       tptr[1].v1,
				       tptr[2].v1,
				       tptr[3].v1);
	      
	const mic_f v2 = gather_4f_zlc(and_mask,
				       tptr[0].v2,
				       tptr[1].v2,
				       tptr[2].v2,
				       tptr[3].v2);

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
	const mic_i rayMask(ray16.mask[rayIndex]);
	const mic_i triMask0( scene->getTriangleMesh( tptr[0].geomID() )->mask );
	const mic_i triMask1( scene->getTriangleMesh( tptr[1].geomID() )->mask );
	const mic_i triMask2( scene->getTriangleMesh( tptr[2].geomID() )->mask );
	const mic_i triMask3( scene->getTriangleMesh( tptr[3].geomID() )->mask );
	const mic_i triMask = gather16i_4i_align(triMask0,triMask1,triMask2,triMask3);
	const mic_m m_ray_mask = (rayMask & triMask) != mic_i::zero();
	m_final &= m_ray_mask;	      
#endif

	//////////////////////////////////////////////////////////////////////////////////////////////////

	/* did the ray hit one of the four triangles? */
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
		    const Triangle1mc  *__restrict__ tri_ptr = tptr + triIndex;
		    const mic_m m_tri = m_dist^(m_dist & (mic_m)((unsigned int)m_dist - 1));
		    const int geomID = tri_ptr->geomID();
		    const int primID = tri_ptr->primID();

		    const mic_f gnormalx(normal[triIndex*4+1]);
		    const mic_f gnormaly(normal[triIndex*4+2]);
		    const mic_f gnormalz(normal[triIndex*4+0]);

                
		    const Geometry* const geom = scene->get(geomID);
		    if (likely(!geom->hasIntersectionFilter16())) 
		      {

			compactustore16f_low(m_tri,&ray16.tfar[rayIndex],min_dist);
			compactustore16f_low(m_tri,&ray16.u[rayIndex],u); 
			compactustore16f_low(m_tri,&ray16.v[rayIndex],v); 
			compactustore16f_low(m_tri,&ray16.Ng.x[rayIndex],gnormalx); 
			compactustore16f_low(m_tri,&ray16.Ng.y[rayIndex],gnormaly); 
			compactustore16f_low(m_tri,&ray16.Ng.z[rayIndex],gnormalz); 
			ray16.geomID[rayIndex] = geomID;
			ray16.primID[rayIndex] = primID;
			max_dist_xyz = min_dist;
			break;
		      }
                
		    if (runIntersectionFilter16(geom,ray16,rayIndex,u,v,min_dist,gnormalx,gnormaly,gnormalz,m_tri,geomID,primID)) {
		      max_dist_xyz = min_dist;
		      break;
		    }
		    m_final ^= m_tri;
		  }
		max_dist_xyz = ray16.tfar[rayIndex];
	      }
	    else
	      {
		STAT3(normal.trav_prim_hits,1,1,1);
		max_dist_xyz  = select(m_final,t,max_dist_xyz);
		const mic_f min_dist = vreduce_min(max_dist_xyz);
		const mic_m m_dist = eq(min_dist,max_dist_xyz);
		const size_t vecIndex = bitscan(toInt(m_dist));
		const size_t triIndex = vecIndex >> 2;

		const Triangle1mc  *__restrict__ tri_ptr = tptr + triIndex;

		const mic_m m_tri = m_dist^(m_dist & (mic_m)((unsigned int)m_dist - 1));

                
		prefetch<PFHINT_L1EX>(&ray16.tfar);  
		prefetch<PFHINT_L1EX>(&ray16.u);
		prefetch<PFHINT_L1EX>(&ray16.v);
		prefetch<PFHINT_L1EX>(&ray16.Ng.x); 
		prefetch<PFHINT_L1EX>(&ray16.Ng.y); 
		prefetch<PFHINT_L1EX>(&ray16.Ng.z); 
		prefetch<PFHINT_L1EX>(&ray16.geomID);
		prefetch<PFHINT_L1EX>(&ray16.primID);

		max_dist_xyz = min_dist;

		const mic_f gnormalx(normal[triIndex*4+1]);
		const mic_f gnormaly(normal[triIndex*4+2]);
		const mic_f gnormalz(normal[triIndex*4+0]);
		  
		compactustore16f_low(m_tri,&ray16.tfar[rayIndex],min_dist);
		compactustore16f_low(m_tri,&ray16.u[rayIndex],u); 
		compactustore16f_low(m_tri,&ray16.v[rayIndex],v); 
		compactustore16f_low(m_tri,&ray16.Ng.x[rayIndex],gnormalx); 
		compactustore16f_low(m_tri,&ray16.Ng.y[rayIndex],gnormaly); 
		compactustore16f_low(m_tri,&ray16.Ng.z[rayIndex],gnormalz); 

		ray16.geomID[rayIndex] = tri_ptr->geomID();
		ray16.primID[rayIndex] = tri_ptr->primID();
	      }
	    return true;
      
	  }
	return false;
      }


      __forceinline static bool occluded1(const size_t rayIndex, 
					  const mic_f &dir_xyz,
					  const mic_f &org_xyz,
					  const mic_f &min_dist_xyz,
					  const mic_f &max_dist_xyz,
					  const mic_i &and_mask,
					  const Ray16& ray16, 
					  mic_m &m_terminated,
					  const Scene     *__restrict__ const scene,
					  const Triangle1mc * __restrict__ const tptr)
      {
	const mic_f zero = mic_f::zero();
	prefetch<PFHINT_L1>(tptr + 0); 
	prefetch<PFHINT_L1>(tptr + 2);

	      
	const mic_f v0 = gather_4f_zlc(and_mask,
				       tptr[0].v0,
				       tptr[1].v0,
				       tptr[2].v0,
				       tptr[3].v0);
	      
	const mic_f v1 = gather_4f_zlc(and_mask,
				       tptr[0].v1,
				       tptr[1].v1,
				       tptr[2].v1,
				       tptr[3].v1);
	      
	const mic_f v2 = gather_4f_zlc(and_mask,
				       tptr[0].v2,
				       tptr[1].v2,
				       tptr[2].v2,
				       tptr[3].v2);

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
	if (unlikely(none(m_aperture))) return false;

	mic_m m_final  = lt(lt(m_aperture,min_dist_xyz,t),t,max_dist_xyz);

#if defined(RTCORE_RAY_MASK)
	const mic_i rayMask(ray16.mask[rayIndex]);
	const mic_i triMask0( scene->getTriangleMesh( tptr[0].geomID() )->mask );
	const mic_i triMask1( scene->getTriangleMesh( tptr[1].geomID() )->mask );
	const mic_i triMask2( scene->getTriangleMesh( tptr[2].geomID() )->mask );
	const mic_i triMask3( scene->getTriangleMesh( tptr[3].geomID() )->mask );
	const mic_i triMask = gather16i_4i_align(triMask0,triMask1,triMask2,triMask3);

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
		const Triangle1mc  *__restrict__ tri_ptr = tptr + triIndex;
		const mic_m m_tri = m_dist^(m_dist & (mic_m)((unsigned int)m_dist - 1));
		const int geomID = tri_ptr->geomID();
		const int primID = tri_ptr->primID();                
		const Geometry* const geom = scene->get(geomID);
		const mic_f gnormalx(normal[triIndex*4+1]);
		const mic_f gnormaly(normal[triIndex*4+2]);
		const mic_f gnormalz(normal[triIndex*4+0]);

		if (likely(!geom->hasOcclusionFilter16())) break;
                
		if (runOcclusionFilter16(geom,(Ray16&)ray16,rayIndex,u,v,min_dist,gnormalx,gnormaly,gnormalz,m_tri,geomID,primID)) 
		  break;

		m_final ^= m_tri; /* clear bit */
	      }
      
	  }
	if (unlikely(any(m_final)))
	  {
	    STAT3(shadow.trav_prim_hits,1,1,1);
	    m_terminated |= mic_m::shift1[rayIndex];
	    return true;
	  }
	return false;
      }

      // ==================================================================================================
      // ==================================================================================================
      // ==================================================================================================


      __forceinline static void intersect16(const mic_m valid_leaf, 
					    const unsigned int items,
					    const mic3f &dir,
					    const mic3f &org,
					    Ray16& ray16, 
					    const Scene     *__restrict__ const scene,
					    const Triangle1mc * __restrict__ tptr)
      {
	const mic_f zero = mic_f::zero();
	const mic_f one  = mic_f::one();
	
	prefetch<PFHINT_L1>((mic_f*)tptr +  0); 
	prefetch<PFHINT_L1>((mic_f*)tptr +  2); 

	for (size_t i=0; i<items; i++,tptr++) 
	  {

	    const Triangle1mc& tri = *tptr;
	    STAT3(normal.trav_prims,1,popcnt(valid_leaf),16);

	    /* load vertices and calculate edges */
	    const mic_f v0 = broadcast4to16f(tri.v0);
	    const mic_f v1 = broadcast4to16f(tri.v1);
	    const mic_f v2 = broadcast4to16f(tri.v2);
	    

	    const mic_f e1 = v0-v1;
	    const mic_f e2 = v2-v0;

	    /* calculate denominator */
	    const mic3f _v0 = mic3f(swizzle<0>(v0),swizzle<1>(v0),swizzle<2>(v0));
	    const mic3f C =  _v0 - org;
	    
	    //const mic3f Ng = mic3f(tri.Ng);
	    const mic_f normal_xyz = lcross_zxy(e1,e2);
	    const mic3f Ng(swBBBB(normal_xyz),swCCCC(normal_xyz),swAAAA(normal_xyz));

	    const mic_f den = dot(Ng,ray16.dir);

	    const mic_f rcp_den = rcp(den);

	    mic_m valid = valid_leaf;

#if defined(RTCORE_BACKFACE_CULLING)
	    
	    valid &= den > zero;
#endif

	    /* perform edge tests */
	    const mic3f R = -cross(C,ray16.dir);
	    const mic3f _e2(swizzle<0>(e2),swizzle<1>(e2),swizzle<2>(e2));
	    const mic_f u = dot(R,_e2)*rcp_den;
	    const mic3f _e1(swizzle<0>(e1),swizzle<1>(e1),swizzle<2>(e1));
	    const mic_f v = dot(R,_e1)*rcp_den;
	    valid = ge(valid,u,zero);
	    valid = ge(valid,v,zero);
	    valid = le(valid,u+v,one);

	    prefetch<PFHINT_L1EX>(&ray16.u);      
	    prefetch<PFHINT_L1EX>(&ray16.v);      
	    prefetch<PFHINT_L1EX>(&ray16.tfar);      


	    if (unlikely(none(valid))) continue;

	    const mic_f dot_C_Ng = dot(C,Ng);
	    const mic_f t = dot_C_Ng * rcp_den;
      
	    /* perform depth test */
	    valid = ge(valid, t,ray16.tnear);
	    valid = ge(valid,ray16.tfar,t);

	    const mic_i geomID = tri.geomID();
	    const mic_i primID = tri.primID();
	    prefetch<PFHINT_L1EX>(&ray16.geomID);      
	    prefetch<PFHINT_L1EX>(&ray16.primID);      
	    prefetch<PFHINT_L1EX>(&ray16.Ng.x);      
	    prefetch<PFHINT_L1EX>(&ray16.Ng.y);      
	    prefetch<PFHINT_L1EX>(&ray16.Ng.z);      

	    /* ray masking test */
#if defined(RTCORE_RAY_MASK)
	    const mic_i tri_mask = scene->getTriangleMesh( tri.geomID() )->mask;
	    valid &= (tri_mask & ray16.mask) != 0;
#endif
	    if (unlikely(none(valid))) continue;

	    /* intersection filter test */
	    if (ENABLE_INTERSECTION_FILTER) 
	      {

		Geometry* geom = ((Scene*)scene)->get(tri.geomID());
		if (unlikely(geom->hasIntersectionFilter16())) {
		  runIntersectionFilter16(valid,geom,ray16,u,v,t,Ng,geomID,primID);
		  continue;
		}
	      }
	    /* update hit information */
	    store16f(valid,(float*)&ray16.u,u);
	    store16f(valid,(float*)&ray16.v,v);
	    store16f(valid,(float*)&ray16.tfar,t);
	    store16i(valid,(float*)&ray16.geomID,geomID);
	    store16i(valid,(float*)&ray16.primID,primID);
	    store16f(valid,(float*)&ray16.Ng.x,Ng.x);
	    store16f(valid,(float*)&ray16.Ng.y,Ng.y);
	    store16f(valid,(float*)&ray16.Ng.z,Ng.z);
	  }

      }

      __forceinline static void occluded16(const mic_m m_valid_leaf, 
					   const unsigned int items,
					   const mic3f &dir,
					   const mic3f &org,
					   Ray16& ray16, 
					   mic_m &m_terminated,
					   const Scene     *__restrict__ const scene,
					   const Triangle1mc * __restrict__ tptr)
      {
	prefetch<PFHINT_L1>((mic_f*)tptr +  0); 
	prefetch<PFHINT_L1>((mic_f*)tptr +  2); 

	const mic_f zero = mic_f::zero();
	mic_m valid_leaf = m_valid_leaf;
	for (size_t i=0; i<items; i++,tptr++) 
	  {
	    const Triangle1mc& tri = *tptr;

	    STAT3(normal.trav_prims,1,popcnt(m_valid_leaf),16);
        
	    /* load vertices and calculate edges */
	    const mic_f v0 = broadcast4to16f(tri.v0);
	    const mic_f v1 = broadcast4to16f(tri.v1);
	    const mic_f v2 = broadcast4to16f(tri.v2);
	    const mic_f e1 = v0-v1;
	    const mic_f e2 = v2-v0;

	    /* calculate denominator */
	    const mic3f _v0 = mic3f(swizzle<0>(v0),swizzle<1>(v0),swizzle<2>(v0));
	    const mic3f C =  _v0 - org;
	    
	    const mic_f normal_xyz = lcross_zxy(e1,e2);
	    const mic3f Ng(swBBBB(normal_xyz),swCCCC(normal_xyz),swAAAA(normal_xyz));

	    const mic_f den = dot(Ng,dir);

	    mic_m valid = valid_leaf;

#if defined(RTCORE_BACKFACE_CULLING)
	    
	    valid &= den > zero;
#endif

	    /* perform edge tests */
	    const mic_f rcp_den = rcp(den);
	    const mic3f R = cross(dir,C);
	    const mic3f _e2(swizzle<0>(e2),swizzle<1>(e2),swizzle<2>(e2));
	    const mic_f u = dot(R,_e2)*rcp_den;
	    const mic3f _e1(swizzle<0>(e1),swizzle<1>(e1),swizzle<2>(e1));
	    const mic_f v = dot(R,_e1)*rcp_den;
	    valid = ge(valid,u,zero);
	    valid = ge(valid,v,zero);
	    valid = le(valid,u+v,one);
	    const mic_f t = dot(C,Ng) * rcp_den;

	    if (unlikely(none(valid))) continue;
      
	    /* perform depth test */
	    valid = ge(valid, t,ray16.tnear);
	    valid = ge(valid,ray16.tfar,t);

	    /* ray masking test */
#if defined(RTCORE_RAY_MASK)
	    const mic_i tri_mask = scene->getTriangleMesh( tri.geomID() )->mask;
	    valid &= (tri_mask & ray16.mask) != 0;
#endif
	    if (unlikely(none(valid))) continue;
	    
	    /* intersection filter test */
	    if (ENABLE_INTERSECTION_FILTER) 
	      {
		const int geomID = tri.geomID();
		const int primID = tri.primID();
		const Geometry* geom = scene->get(geomID);
		if (unlikely(geom->hasOcclusionFilter16()))
		  valid = runOcclusionFilter16(valid,geom,(Ray16&)ray16,u,v,t,Ng,geomID,primID);
	      }
	    /* update occlusion */
	    m_terminated |= valid;
	  valid_leaf &= ~valid;
	  if (unlikely(none(valid_leaf))) break;
	}
    }

  };
}
