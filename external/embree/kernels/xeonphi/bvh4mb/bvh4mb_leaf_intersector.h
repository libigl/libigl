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

#include "bvh4mb.h"

namespace embree
{
  static __aligned(64) int zlc4[4] = {0xffffffff,0xffffffff,0xffffffff,0};

  struct Triangle1mbLeafIntersector
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
      const mic_f time     = broadcast1to16f(&ray.time);
      const mic_f one_time = (mic_f::one() - time);

      const BVH4mb::Triangle01* tptr  = (BVH4mb::Triangle01*) curNode.leaf<8>(accel);
      prefetch<PFHINT_L2>((mic_f*)tptr +  0); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  1); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  2); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  3); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  4); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  5); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  6); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  7); 

      const mic_i and_mask = broadcast4to16i(zlc4);
	      
      const mic_f v0_t0 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t0.v0,
					(float*)&tptr[1].t0.v0,
					(float*)&tptr[2].t0.v0,
					(float*)&tptr[3].t0.v0);
	      
      const mic_f v1_t0 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t0.v1,
					(float*)&tptr[1].t0.v1,
					(float*)&tptr[2].t0.v1,
					(float*)&tptr[3].t0.v1);
	      
      const mic_f v2_t0 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t0.v2,
					(float*)&tptr[1].t0.v2,
					(float*)&tptr[2].t0.v2,
					(float*)&tptr[3].t0.v2);

      const mic_f v0_t1 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t1.v0,
					(float*)&tptr[1].t1.v0,
					(float*)&tptr[2].t1.v0,
					(float*)&tptr[3].t1.v0);
	      
      const mic_f v1_t1 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t1.v1,
					(float*)&tptr[1].t1.v1,
					(float*)&tptr[2].t1.v1,
					(float*)&tptr[3].t1.v1);
	      
      const mic_f v2_t1 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t1.v2,
					(float*)&tptr[1].t1.v2,
					(float*)&tptr[2].t1.v2,
					(float*)&tptr[3].t1.v2);

      const mic_f v0 = v0_t0 * one_time + time * v0_t1;
      const mic_f v1 = v1_t0 * one_time + time * v1_t1;
      const mic_f v2 = v2_t0 * one_time + time * v2_t1;

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

      max_dist_xyz  = select(m_final,t,max_dist_xyz);
		    
      //////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(RTCORE_RAY_MASK)
      const mic_i rayMask(ray.mask);
      const mic_i triMask = getTriMasks(tptr); 
      const mic_m m_ray_mask = (rayMask & triMask) != mic_i::zero();
      m_final &= m_ray_mask;	      
#endif


      /* did the ray hot one of the four triangles? */
      if (unlikely(any(m_final)))
	{
	  const mic_f min_dist = vreduce_min(max_dist_xyz);
	  const mic_m m_dist = eq(min_dist,max_dist_xyz);

	  prefetch<PFHINT_L1EX>((mic_f*)&ray + 0);
	  prefetch<PFHINT_L1EX>((mic_f*)&ray + 1);

	  const size_t vecIndex = bitscan(toInt(m_dist));
	  const size_t triIndex = vecIndex >> 2;

	  const BVH4mb::Triangle01  *__restrict__ tri_ptr = tptr + triIndex;

	  const mic_m m_tri = m_dist^(m_dist & (mic_m)((unsigned int)m_dist - 1));

	  const mic_f gnormalz = swAAAA(normal);
	  const mic_f gnormalx = swBBBB(normal);
	  const mic_f gnormaly = swCCCC(normal);
		  
	  max_dist_xyz = min_dist;

	  compactustore16f_low(m_tri,&ray.tfar,min_dist);
	  compactustore16f_low(m_tri,&ray.u,u); 
	  compactustore16f_low(m_tri,&ray.v,v); 
	  compactustore16f_low(m_tri,&ray.Ng.x,gnormalx); 
	  compactustore16f_low(m_tri,&ray.Ng.y,gnormaly); 
	  compactustore16f_low(m_tri,&ray.Ng.z,gnormalz); 

	  ray.geomID = tri_ptr->t0.geomID();
	  ray.primID = tri_ptr->t0.primID();
	  return true;
	}
      return false;
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
      const mic_f time     = broadcast1to16f(&ray.time);
      const mic_f one_time = (mic_f::one() - time);
      const BVH4mb::Triangle01* tptr  = (BVH4mb::Triangle01*) curNode.leaf<8>(accel);

      prefetch<PFHINT_L2>((mic_f*)tptr +  0); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  1); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  2); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  3); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  4); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  5); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  6); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  7); 

      const mic_i and_mask = broadcast4to16i(zlc4);
	      
      const mic_f v0_t0 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t0.v0,
					(float*)&tptr[1].t0.v0,
					(float*)&tptr[2].t0.v0,
					(float*)&tptr[3].t0.v0);
	      
      const mic_f v1_t0 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t0.v1,
					(float*)&tptr[1].t0.v1,
					(float*)&tptr[2].t0.v1,
					(float*)&tptr[3].t0.v1);
	      
      const mic_f v2_t0 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t0.v2,
					(float*)&tptr[1].t0.v2,
					(float*)&tptr[2].t0.v2,
					(float*)&tptr[3].t0.v2);

      const mic_f v0_t1 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t1.v0,
					(float*)&tptr[1].t1.v0,
					(float*)&tptr[2].t1.v0,
					(float*)&tptr[3].t1.v0);
	      
      const mic_f v1_t1 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t1.v1,
					(float*)&tptr[1].t1.v1,
					(float*)&tptr[2].t1.v1,
					(float*)&tptr[3].t1.v1);
	      
      const mic_f v2_t1 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t1.v2,
					(float*)&tptr[1].t1.v2,
					(float*)&tptr[2].t1.v2,
					(float*)&tptr[3].t1.v2);


      const mic_f v0 = v0_t0 * one_time + time * v0_t1;
      const mic_f v1 = v1_t0 * one_time + time * v1_t1;
      const mic_f v2 = v2_t0 * one_time + time * v2_t1;

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
      const mic_f t = rcp_den*nom;

      if (unlikely(none(m_aperture))) return false;

      mic_m m_final  = lt(lt(m_aperture,min_dist_xyz,t),t,max_dist_xyz);

#if defined(RTCORE_RAY_MASK)
      const mic_i rayMask(ray.mask);
      const mic_i triMask = getTriMasks(tptr); 
      const mic_m m_ray_mask = (rayMask & triMask) != mic_i::zero();
      m_final &= m_ray_mask;	      
#endif

      return any(m_final);
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
      const mic_f time     = broadcast1to16f(&ray16.time[rayIndex]);
      const mic_f one_time = (mic_f::one() - time);

      const BVH4mb::Triangle01* tptr  = (BVH4mb::Triangle01*) curNode.leaf<8>(accel);
	      
      prefetch<PFHINT_L1>((mic_f*)tptr +  0); 
      prefetch<PFHINT_L1>((mic_f*)tptr +  1); 
      prefetch<PFHINT_L1>((mic_f*)tptr +  2); 
      prefetch<PFHINT_L1>((mic_f*)tptr +  3); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  4); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  5); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  6); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  7); 

      const mic_i and_mask = broadcast4to16i(zlc4);
	     

      const mic_f v0_t0 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t0.v0,
					(float*)&tptr[1].t0.v0,
					(float*)&tptr[2].t0.v0,
					(float*)&tptr[3].t0.v0);
	      
      const mic_f v1_t0 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t0.v1,
					(float*)&tptr[1].t0.v1,
					(float*)&tptr[2].t0.v1,
					(float*)&tptr[3].t0.v1);
	      
      const mic_f v2_t0 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t0.v2,
					(float*)&tptr[1].t0.v2,
					(float*)&tptr[2].t0.v2,
					(float*)&tptr[3].t0.v2);

      const mic_f v0_t1 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t1.v0,
					(float*)&tptr[1].t1.v0,
					(float*)&tptr[2].t1.v0,
					(float*)&tptr[3].t1.v0);
	      
      const mic_f v1_t1 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t1.v1,
					(float*)&tptr[1].t1.v1,
					(float*)&tptr[2].t1.v1,
					(float*)&tptr[3].t1.v1);
	      
      const mic_f v2_t1 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t1.v2,
					(float*)&tptr[1].t1.v2,
					(float*)&tptr[2].t1.v2,
					(float*)&tptr[3].t1.v2);

      const mic_f v0 = v0_t0 * one_time + time * v0_t1;
      const mic_f v1 = v1_t0 * one_time + time * v1_t1;
      const mic_f v2 = v2_t0 * one_time + time * v2_t1;

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

      max_dist_xyz  = select(m_final,t,max_dist_xyz);
		    
      //////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(RTCORE_RAY_MASK)
      const mic_i rayMask(ray16.mask[rayIndex]);
      const mic_i triMask = getTriMasks(tptr); 
      const mic_m m_ray_mask = (rayMask & triMask) != mic_i::zero();
      m_final &= m_ray_mask;	      
#endif


      /* did the ray hot one of the four triangles? */
      if (unlikely(any(m_final)))
	{
	  const mic_f min_dist = vreduce_min(max_dist_xyz);
	  const mic_m m_dist = eq(min_dist,max_dist_xyz);

	  const size_t vecIndex = bitscan(toInt(m_dist));
	  const size_t triIndex = vecIndex >> 2;

	  const BVH4mb::Triangle01  *__restrict__ tri_ptr = tptr + triIndex;

	  const mic_m m_tri = m_dist^(m_dist & (mic_m)((unsigned int)m_dist - 1));
		  
	  const mic_f gnormalz = swAAAA(normal);
	  const mic_f gnormalx = swBBBB(normal);
	  const mic_f gnormaly = swCCCC(normal);

	  prefetch<PFHINT_L1EX>(&ray16.tfar);  
	  prefetch<PFHINT_L1EX>(&ray16.u);
	  prefetch<PFHINT_L1EX>(&ray16.v);
	  prefetch<PFHINT_L1EX>(&ray16.Ng.x); 
	  prefetch<PFHINT_L1EX>(&ray16.Ng.y); 
	  prefetch<PFHINT_L1EX>(&ray16.Ng.z); 
	  prefetch<PFHINT_L1EX>(&ray16.geomID);
	  prefetch<PFHINT_L1EX>(&ray16.primID);

	  max_dist_xyz = min_dist;
		  
	  compactustore16f_low(m_tri,&ray16.tfar[rayIndex],min_dist);
	  compactustore16f_low(m_tri,&ray16.u[rayIndex],u); 
	  compactustore16f_low(m_tri,&ray16.v[rayIndex],v); 
	  compactustore16f_low(m_tri,&ray16.Ng.x[rayIndex],gnormalx); 
	  compactustore16f_low(m_tri,&ray16.Ng.y[rayIndex],gnormaly); 
	  compactustore16f_low(m_tri,&ray16.Ng.z[rayIndex],gnormalz); 

	  ray16.geomID[rayIndex] = tri_ptr->t0.geomID();
	  ray16.primID[rayIndex] = tri_ptr->t0.primID();
	  return true;
	}
      return false;
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
      const mic_f time     = broadcast1to16f(&ray16.time[rayIndex]);
      const mic_f one_time = (mic_f::one() - time);

      const BVH4mb::Triangle01* tptr  = (BVH4mb::Triangle01*) curNode.leaf<8>(accel);

      prefetch<PFHINT_L1>((mic_f*)tptr +  0); 
      prefetch<PFHINT_L1>((mic_f*)tptr +  1); 
      prefetch<PFHINT_L1>((mic_f*)tptr +  2); 
      prefetch<PFHINT_L1>((mic_f*)tptr +  3); 

      const mic_i and_mask = broadcast4to16i(zlc4);
	      
      const mic_f v0_t0 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t0.v0,
					(float*)&tptr[1].t0.v0,
					(float*)&tptr[2].t0.v0,
					(float*)&tptr[3].t0.v0);
	      
      const mic_f v1_t0 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t0.v1,
					(float*)&tptr[1].t0.v1,
					(float*)&tptr[2].t0.v1,
					(float*)&tptr[3].t0.v1);
	      
      const mic_f v2_t0 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t0.v2,
					(float*)&tptr[1].t0.v2,
					(float*)&tptr[2].t0.v2,
					(float*)&tptr[3].t0.v2);


      prefetch<PFHINT_L2>((mic_f*)tptr +  4); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  5); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  6); 
      prefetch<PFHINT_L2>((mic_f*)tptr +  7); 

      const mic_f v0_t1 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t1.v0,
					(float*)&tptr[1].t1.v0,
					(float*)&tptr[2].t1.v0,
					(float*)&tptr[3].t1.v0);
	      
      const mic_f v1_t1 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t1.v1,
					(float*)&tptr[1].t1.v1,
					(float*)&tptr[2].t1.v1,
					(float*)&tptr[3].t1.v1);
	      
      const mic_f v2_t1 = gather_4f_zlc(and_mask,
					(float*)&tptr[0].t1.v2,
					(float*)&tptr[1].t1.v2,
					(float*)&tptr[2].t1.v2,
					(float*)&tptr[3].t1.v2);

      const mic_f v0 = v0_t0 * one_time + time * v0_t1;
      const mic_f v1 = v1_t0 * one_time + time * v1_t1;
      const mic_f v2 = v2_t0 * one_time + time * v2_t1;

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
      const mic_i triMask = getTriMasks(tptr); 
      const mic_m m_ray_mask = (rayMask & triMask) != mic_i::zero();
      m_final &= m_ray_mask;	      
#endif

      if (unlikely(any(m_final)))
	{
	  m_terminated |= mic_m::shift1[rayIndex];
	  return true;
	}
      return false;
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

	const mic_f time     = ray16.time;
	const mic_f one_time = (mic_f::one() - time);

	unsigned int items; 
	const BVH4mb::Triangle01* tris  = (BVH4mb::Triangle01*) curNode.leaf<8>(accel,items);

	const mic_f zero = mic_f::zero();
	const mic_f one  = mic_f::one();

	prefetch<PFHINT_L1>((mic_f*)tris +  0); 
	prefetch<PFHINT_L2>((mic_f*)tris +  1); 
	prefetch<PFHINT_L2>((mic_f*)tris +  2); 
	prefetch<PFHINT_L2>((mic_f*)tris +  3); 
	prefetch<PFHINT_L2>((mic_f*)tris +  4); 
	prefetch<PFHINT_L2>((mic_f*)tris +  5); 
	prefetch<PFHINT_L2>((mic_f*)tris +  6); 
	prefetch<PFHINT_L2>((mic_f*)tris +  7); 

	for (size_t i=0; i<items; i++) 
	  {
	    const Triangle1& tri_t0 = tris[i].t0;
	    const Triangle1& tri_t1 = tris[i].t1;

	    prefetch<PFHINT_L1>(&tris[i+1].t0); 
	    prefetch<PFHINT_L1>(&tris[i+1].t1); 

	    STAT3(normal.trav_prims,1,popcnt(m_valid_leaf),16);
        
	    /* load vertices and calculate edges */
	    const mic3f v0_t0( broadcast1to16f(&tri_t0.v0.x), broadcast1to16f(&tri_t0.v0.y), broadcast1to16f(&tri_t0.v0.z) );
	    const mic3f v0_t1( broadcast1to16f(&tri_t1.v0.x), broadcast1to16f(&tri_t1.v0.y), broadcast1to16f(&tri_t1.v0.z) );
	    const mic3f v0 = v0_t0 * one_time + time * v0_t1;
	    const mic3f v1_t0( broadcast1to16f(&tri_t0.v1.x), broadcast1to16f(&tri_t0.v1.y), broadcast1to16f(&tri_t0.v1.z) );
	    const mic3f v1_t1( broadcast1to16f(&tri_t1.v1.x), broadcast1to16f(&tri_t1.v1.y), broadcast1to16f(&tri_t1.v1.z) );
	    const mic3f v1 = v1_t0 * one_time + time * v1_t1;
	    const mic3f v2_t0( broadcast1to16f(&tri_t0.v2.x), broadcast1to16f(&tri_t0.v2.y), broadcast1to16f(&tri_t0.v2.z) );
	    const mic3f v2_t1( broadcast1to16f(&tri_t1.v2.x), broadcast1to16f(&tri_t1.v2.y), broadcast1to16f(&tri_t1.v2.z) );
	    const mic3f v2 = v2_t0 * one_time + time * v2_t1;

	    const mic3f e1 = v0-v1;
	    const mic3f e2 = v2-v0;

	    const mic3f Ng = cross(e1,e2);

	    /* calculate denominator */
	    const mic3f C =  v0 - org;
	    
	    const mic_f den = dot(Ng,dir);

	    mic_m valid = m_valid_leaf;

#if defined(RTCORE_BACKFACE_CULLING)
	    
	    valid &= den > zero;
#endif

	    /* perform edge tests */
	    const mic_f rcp_den = rcp(den);
	    const mic3f R = cross(dir,C);
	    const mic_f u = dot(R,e2)*rcp_den;
	    const mic_f v = dot(R,e1)*rcp_den;
	    valid = ge(valid,u,zero);
	    valid = ge(valid,v,zero);
	    valid = le(valid,u+v,one);
	    prefetch<PFHINT_L1EX>(&ray16.u);      
	    prefetch<PFHINT_L1EX>(&ray16.v);      
	    prefetch<PFHINT_L1EX>(&ray16.tfar);      
	    const mic_f t = dot(C,Ng) * rcp_den;

	    if (unlikely(none(valid))) continue;
      
	    /* perform depth test */
	    valid = ge(valid, t,ray16.tnear);
	    valid = ge(valid,ray16.tfar,t);

	    const mic_i geomID = tri_t0.geomID();
	    const mic_i primID = tri_t0.primID();
	    prefetch<PFHINT_L1EX>(&ray16.geomID);      
	    prefetch<PFHINT_L1EX>(&ray16.primID);      
	    prefetch<PFHINT_L1EX>(&ray16.Ng.x);      
	    prefetch<PFHINT_L1EX>(&ray16.Ng.y);      
	    prefetch<PFHINT_L1EX>(&ray16.Ng.z);      

	    /* ray masking test */
#if defined(RTCORE_RAY_MASK)
	    valid &= (mic_i(tri_t0.mask()) & ray16.mask) != 0;
#endif
	    if (unlikely(none(valid))) continue;
        
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

      __forceinline static void occluded16(BVH4i::NodeRef curNode,
					   const mic_m m_valid_leaf_active, 
					   const mic3f &dir,
					   const mic3f &org,
					   Ray16& ray16, 
					   mic_m &m_terminated,					    
					   const void *__restrict__ const accel,
					   const Scene     *__restrict__ const geometry)
      {
	mic_m m_valid_leaf = m_valid_leaf_active;

	const mic_f time     = ray16.time;
	const mic_f one_time = (mic_f::one() - time);

	unsigned int items; 
	const BVH4mb::Triangle01* tris  = (BVH4mb::Triangle01*) curNode.leaf<8>(accel,items);

	prefetch<PFHINT_L1>((mic_f*)tris +  0); 
	prefetch<PFHINT_L2>((mic_f*)tris +  1); 
	prefetch<PFHINT_L2>((mic_f*)tris +  2); 
	prefetch<PFHINT_L2>((mic_f*)tris +  3); 
	prefetch<PFHINT_L2>((mic_f*)tris +  4); 
	prefetch<PFHINT_L2>((mic_f*)tris +  5); 
	prefetch<PFHINT_L2>((mic_f*)tris +  6); 
	prefetch<PFHINT_L2>((mic_f*)tris +  7); 

	const mic_f zero = mic_f::zero();

	for (size_t i=0; i<items; i++) 
	  {
	    const Triangle1& tri_t0 = tris[i].t0;
	    const Triangle1& tri_t1 = tris[i].t1;

	    prefetch<PFHINT_L1>(&tris[i+1].t0); 
	    prefetch<PFHINT_L1>(&tris[i+1].t1); 

	    STAT3(normal.trav_prims,1,popcnt(m_valid_leaf_active),16);
        
	    /* load vertices and calculate edges */
	    const mic3f v0_t0( broadcast1to16f(&tri_t0.v0.x), broadcast1to16f(&tri_t0.v0.y), broadcast1to16f(&tri_t0.v0.z) );
	    const mic3f v0_t1( broadcast1to16f(&tri_t1.v0.x), broadcast1to16f(&tri_t1.v0.y), broadcast1to16f(&tri_t1.v0.z) );
	    const mic3f v0 = v0_t0 * one_time + time * v0_t1;
	    const mic3f v1_t0( broadcast1to16f(&tri_t0.v1.x), broadcast1to16f(&tri_t0.v1.y), broadcast1to16f(&tri_t0.v1.z) );
	    const mic3f v1_t1( broadcast1to16f(&tri_t1.v1.x), broadcast1to16f(&tri_t1.v1.y), broadcast1to16f(&tri_t1.v1.z) );
	    const mic3f v1 = v1_t0 * one_time + time * v1_t1;
	    const mic3f v2_t0( broadcast1to16f(&tri_t0.v2.x), broadcast1to16f(&tri_t0.v2.y), broadcast1to16f(&tri_t0.v2.z) );
	    const mic3f v2_t1( broadcast1to16f(&tri_t1.v2.x), broadcast1to16f(&tri_t1.v2.y), broadcast1to16f(&tri_t1.v2.z) );
	    const mic3f v2 = v2_t0 * one_time + time * v2_t1;

	    const mic3f e1 = v0-v1;
	    const mic3f e2 = v2-v0;

	    const mic3f Ng = cross(e1,e2);

	    /* calculate denominator */
	    const mic3f C =  v0 - org;
	    
	    const mic_f den = dot(Ng,dir);

	    mic_m valid = m_valid_leaf;

#if defined(RTCORE_BACKFACE_CULLING)
	    
	    valid &= den > zero;
#endif

	    /* perform edge tests */
	    const mic_f rcp_den = rcp(den);
	    const mic3f R = cross(dir,C);
	    const mic_f u = dot(R,e2)*rcp_den;
	    const mic_f v = dot(R,e1)*rcp_den;
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
	    valid &= (mic_i(tri_t0.mask()) & ray16.mask) != 0;
#endif
	    if (unlikely(none(valid))) continue;
	    
	    /* update occlusion */
	    m_terminated |= valid;
	    m_valid_leaf &= ~valid;
	    if (unlikely(none(m_valid_leaf))) break;
	  }

      }

    };

};
