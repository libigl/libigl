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

#include "bezier1i.h"
#include "common/ray16.h"
#include "geometry/filter.h"

using namespace std;

namespace embree
{
  typedef LinearSpace3<mic3f> LinearSpace_mic3f;
    
  /*! Intersector for a single ray from a ray packet with a bezier curve. */

  struct __aligned(64) Precalculations 
  {
     /* __forceinline Precalculations (const Ray16& ray, const size_t k)  */
     /*   : ray_space(frame(Vec3fa(ray.dir.x[k],ray.dir.y[k],ray.dir.z[k])).transposed()) {} // FIXME: works only with normalized ray direction  */

    __forceinline Precalculations (const LinearSpace_mic3f& ls16, const mic_f &rcp_length, const size_t k)
      : ray_space(ls16.vx.x[k],ls16.vy.x[k],ls16.vz.x[k],
		  ls16.vx.y[k],ls16.vy.y[k],ls16.vz.y[k],
		  ls16.vx.z[k],ls16.vy.z[k],ls16.vz.z[k]),
      inv_ray_length(rcp_length[k])
      {}
    __aligned(64) LinearSpace3fa ray_space;
    mic_f inv_ray_length;
  };

  template< bool ENABLE_INTERSECTION_FILTER>
  struct Bezier1iIntersector16
  {
    typedef Bezier1i Primitive;


    static __forceinline mic4f eval16(const mic_f &p0123,
				      const mic_f &c0,
				      const mic_f &c1,
				      const mic_f &c2,
				      const mic_f &c3)
    {
#if 1
      const mic_f p0 = permute<0>(p0123);
      const mic_f p1 = permute<1>(p0123);
      const mic_f p2 = permute<2>(p0123);
      const mic_f p3 = permute<3>(p0123);

      const mic_f x = c0 * swAAAA(p0) + c1 * swAAAA(p1) + c2 * swAAAA(p2) + c3 * swAAAA(p3);
      const mic_f y = c0 * swBBBB(p0) + c1 * swBBBB(p1) + c2 * swBBBB(p2) + c3 * swBBBB(p3);
      const mic_f z = c0 * swCCCC(p0) + c1 * swCCCC(p1) + c2 * swCCCC(p2) + c3 * swCCCC(p3);
      const mic_f w = c0 * swDDDD(p0) + c1 * swDDDD(p1) + c2 * swDDDD(p2) + c3 * swDDDD(p3);
#else
      const mic_f x = madd(c0,mic_f(p0123[0]),madd(c1,mic_f(p0123[4]),madd(c2,mic_f(p0123[8]),c3* mic_f(p0123[12]))));
      const mic_f y = madd(c0,mic_f(p0123[1]),madd(c1,mic_f(p0123[5]),madd(c2,mic_f(p0123[9]),c3* mic_f(p0123[13]))));
      const mic_f z = madd(c0,mic_f(p0123[2]),madd(c1,mic_f(p0123[6]),madd(c2,mic_f(p0123[10]),c3* mic_f(p0123[14]))));
      const mic_f w = madd(c0,mic_f(p0123[3]),madd(c1,mic_f(p0123[7]),madd(c2,mic_f(p0123[11]),c3* mic_f(p0123[15]))));
#endif
      return mic4f(x,y,z,w);
    }

    static __forceinline void eval(const float t, const mic_f &p0123, mic_f& point, mic_f& tangent)
    {
      const mic_f t0 = mic_f(1.0f) - mic_f(t), t1 = mic_f(t);

      const Vec3fa *__restrict__ const p = (Vec3fa*)&p0123;

      const mic_f p00 = broadcast4to16f((float*)&p[0]);
      const mic_f p01 = broadcast4to16f((float*)&p[1]);
      const mic_f p02 = broadcast4to16f((float*)&p[2]);
      const mic_f p03 = broadcast4to16f((float*)&p[3]);

      const mic_f p10 = p00 * t0 + p01 * t1;
      const mic_f p11 = p01 * t0 + p02 * t1;
      const mic_f p12 = p02 * t0 + p03 * t1;
      const mic_f p20 = p10 * t0 + p11 * t1;
      const mic_f p21 = p11 * t0 + p12 * t1;
      const mic_f p30 = p20 * t0 + p21 * t1;

      point = p30;
      tangent = p21-p20;
    }

    static __forceinline bool intersect(const mic_f& pre_vx, 
					const mic_f& pre_vy, 
					const mic_f& pre_vz, 					
					const mic_f& inv_ray_length,
					Ray16& ray, 
					const mic_f &dir_xyz,
					const mic_f &org_xyz,
					const size_t k, 
					const Bezier1i& curve_in, 
					const void* geom)
    {
      STAT3(normal.trav_prims,1,1,1);

      const mic_f c0 = load16f(&coeff01[0]);
      const mic_f c1 = load16f(&coeff01[1]);
      const mic_f c2 = load16f(&coeff01[2]);
      const mic_f c3 = load16f(&coeff01[3]);

      const mic_f zero = mic_f::zero();
      const mic_f one  = mic_f::one();

      prefetch<PFHINT_L1>(curve_in.p + 0);
      prefetch<PFHINT_L1>(curve_in.p + 3);


      const mic_f p0123 = uload16f((float*)curve_in.p);

      const mic_f p0123_org = p0123 - org_xyz;

      const mic_f p0123_2D = select(0x7777,pre_vx * swAAAA(p0123_org) + pre_vy * swBBBB(p0123_org) + pre_vz * swCCCC(p0123_org),p0123);



      const mic4f p0 = eval16(p0123_2D,c0,c1,c2,c3);
      
      const mic_f last_x = mic_f(p0123_2D[12 + 0]);
      const mic_f last_y = mic_f(p0123_2D[13 + 0]);
      const mic_f last_z = mic_f(p0123_2D[14 + 0]);
      const mic_f last_w = mic_f(p0123_2D[15 + 0]);

      const mic4f p1(align_shift_right<1>(last_x,p0[0]),  
       		     align_shift_right<1>(last_y,p0[1]), 
      		     align_shift_right<1>(last_z,p0[2]),  
       		     align_shift_right<1>(last_w,p0[3]));


      /* approximative intersection with cone */
      const mic4f v = p1-p0;
      const mic4f w = -p0;
      const mic_f d0 = w.x*v.x + w.y*v.y;
      const mic_f d1 = v.x*v.x + v.y*v.y;
      const mic_f u = clamp(d0*rcp_nr(d1),zero,one);
      const mic4f p = p0 + u*v;
      const mic_f t = p.z * inv_ray_length;
      const mic_f d2 = p.x*p.x + p.y*p.y; 
      const mic_f r = p.w;
      const mic_f r2 = r*r;

      mic_m valid = le(d2,r2);
      valid = lt(valid,mic_f(ray.tnear[k]),t);
      valid = lt(valid,t,mic_f(ray.tfar[k]));

      if (unlikely(none(valid))) return false;
      STAT3(normal.trav_prim_hits,1,1,1);

      /* ray masking test */
#if defined(RTCORE_RAY_MASK)
      BezierCurves* g = ((Scene*)geom)->getBezierCurves(curve_in.geomID);
      if (unlikely(g->mask & ray.mask[k]) == 0) return false;
#endif  

      /* update hit information */
      if (ENABLE_INTERSECTION_FILTER) 
	{
	  const Geometry* const gg = ((Scene*)geom)->get(curve_in.geomID);
	  if (unlikely(gg->hasIntersectionFilter16())) 
	    {
	      while(any(valid)) 
		{
		  const float one_over_width = 1.0f/16.0f;
		  unsigned int i = select_min(valid,t);
		  float uu = (float(i)+u[i])*one_over_width; 
		  mic_f P,T;
		  eval(uu,p0123,P,T);
		  assert( T != mic_f::zero() );

		  if (runIntersectionFilter16(gg,ray,k,mic_f(uu),mic_f(0.0f),mic_f(t[i]),mic_f(T[0]),mic_f(T[1]),mic_f(T[2]),(mic_m)((unsigned int)1 << i),curve_in.geomID,curve_in.primID))
		    break;
		  valid ^= (unsigned int)1 << i;
		}
	      if (unlikely(none(valid))) return false;
	    }
	}

      const float one_over_width = 1.0f/16.0f;
      unsigned int i = select_min(valid,t);
      float uu = (float(i)+u[i])*one_over_width; 
      mic_f P,T;
      eval(uu,p0123,P,T);
      assert( T != mic_f::zero() );
      assert( isfinite(T[0]) );
      assert( isfinite(T[1]) );
      assert( isfinite(T[2]) );

      ray.update(1,k,mic_f(t[i]),mic_f(uu),mic_f::zero(),swAAAA(T),swBBBB(T),swCCCC(T),curve_in.geomID,curve_in.primID);

      return true;
    }

    static __forceinline bool occluded(const mic_f& pre_vx, 
				       const mic_f& pre_vy, 
				       const mic_f& pre_vz, 	
				       const mic_f& inv_ray_length,				
				       const Ray16& ray, 
				       const mic_f &dir_xyz,
				       const mic_f &org_xyz,
				       const size_t k, 
				       const Bezier1i& curve_in, 
				       const void* geom) 
    {
      STAT3(shadow.trav_prims,1,1,1);
      const mic_f zero = mic_f::zero();
      const mic_f one  = mic_f::one();

      prefetch<PFHINT_L1>(curve_in.p + 0);
      prefetch<PFHINT_L1>(curve_in.p + 3);

      const mic_f p0123 = uload16f((float*)curve_in.p);
      const mic_f p0123_org = p0123 - org_xyz;

      const mic_f p0123_2D = select(0x7777,pre_vx * swAAAA(p0123_org) + pre_vy * swBBBB(p0123_org) + pre_vz * swCCCC(p0123_org),p0123);


      const mic_f c0 = load16f(&coeff01[0]);
      const mic_f c1 = load16f(&coeff01[1]);
      const mic_f c2 = load16f(&coeff01[2]);
      const mic_f c3 = load16f(&coeff01[3]);

      const mic4f p0 = eval16(p0123_2D,c0,c1,c2,c3);

      const mic_f last_x = mic_f(p0123_2D[12 + 0]);
      const mic_f last_y = mic_f(p0123_2D[13 + 0]);
      const mic_f last_z = mic_f(p0123_2D[14 + 0]);
      const mic_f last_w = mic_f(p0123_2D[15 + 0]);

      const mic4f p1(align_shift_right<1>(last_x,p0[0]),  
       		     align_shift_right<1>(last_y,p0[1]), 
      		     align_shift_right<1>(last_z,p0[2]),  
       		     align_shift_right<1>(last_w,p0[3]));


      const float one_over_width = 1.0f/16.0f;


      /* approximative intersection with cone */
      const mic4f v = p1-p0;
      const mic4f w = -p0;
      const mic_f d0 = w.x*v.x + w.y*v.y;
      const mic_f d1 = v.x*v.x + v.y*v.y;
      const mic_f u = clamp(d0*rcp(d1),zero,one);
      const mic4f p = p0 + u*v;
      const mic_f t = p.z * inv_ray_length;
      const mic_f d2 = p.x*p.x + p.y*p.y; 
      const mic_f r = p.w;
      const mic_f r2 = r*r;
      mic_m valid = le(d2,r2);
      valid = lt(valid,mic_f(ray.tnear[k]),t);
      valid = lt(valid,t,mic_f(ray.tfar[k]));


      if (unlikely(none(valid))) return false;

      STAT3(shadow.trav_prim_hits,1,1,1);

      /* ray masking test */
#if defined(RTCORE_RAY_MASK)
      BezierCurves* g = ((Scene*)geom)->getBezierCurves(curve_in.geomID);
      if (unlikely(g->mask & ray.mask[k]) == 0) return false;
#endif  


      if (ENABLE_INTERSECTION_FILTER) 
	{
	  const Geometry* const gg = ((Scene*)geom)->get(curve_in.geomID);
	  if (likely(gg->hasOcclusionFilter16())) 
	    {
	      while(any(valid)) 
		{
		  unsigned int i = select_min(valid,t);
		  float uu = (float(i)+u[i])*one_over_width; 
		  mic_f P,T;
		  eval(uu,p0123,P,T);
		  assert( T != mic_f::zero() );

		  if (runOcclusionFilter16(gg,(Ray16&)ray,k,mic_f(uu),mic_f(0.0f),mic_f(t[i]),mic_f(T[0]),mic_f(T[1]),mic_f(T[2]),(mic_m)((unsigned int)1 << i),curve_in.geomID,curve_in.primID))
		    return true;
		  valid ^= (unsigned int)1 << i;
		}
	    }
	  return false;
	}
      return true;
    }


    // ==================================================================
    // ==================================================================
    // ==================================================================


    static __forceinline bool intersect(const mic_f& pre_vx, 
					const mic_f& pre_vy, 
					const mic_f& pre_vz, 					
					const mic_f& inv_ray_length,
					Ray& ray, 
					const mic_f &dir_xyz,
					const mic_f &org_xyz,
					const Bezier1i& curve_in, 
					const void* geom)
    {
      STAT3(normal.trav_prims,1,1,1);

      const mic_f zero = mic_f::zero();
      const mic_f one  = mic_f::one();

      prefetch<PFHINT_L1>(curve_in.p + 0);
      prefetch<PFHINT_L1>(curve_in.p + 3);

      const mic_f p0123 = uload16f((float*)curve_in.p);
      const mic_f p0123_org = p0123 - org_xyz;

      const mic_f p0123_2D = select(0x7777,pre_vx * swAAAA(p0123_org) + pre_vy * swBBBB(p0123_org) + pre_vz * swCCCC(p0123_org),p0123);


      const mic_f c0 = load16f(&coeff01[0]);
      const mic_f c1 = load16f(&coeff01[1]);
      const mic_f c2 = load16f(&coeff01[2]);
      const mic_f c3 = load16f(&coeff01[3]);

      const mic4f p0 = eval16(p0123_2D,c0,c1,c2,c3);

      const mic_f last_x = mic_f(p0123_2D[12 + 0]);
      const mic_f last_y = mic_f(p0123_2D[13 + 0]);
      const mic_f last_z = mic_f(p0123_2D[14 + 0]);
      const mic_f last_w = mic_f(p0123_2D[15 + 0]);

      const mic4f p1(align_shift_right<1>(last_x,p0[0]),  
       		     align_shift_right<1>(last_y,p0[1]), 
      		     align_shift_right<1>(last_z,p0[2]),  
       		     align_shift_right<1>(last_w,p0[3]));

      const float one_over_width = 1.0f/16.0f;


      /* approximative intersection with cone */
      const mic4f v = p1-p0;
      const mic4f w = -p0;
      const mic_f d0 = w.x*v.x + w.y*v.y;
      const mic_f d1 = v.x*v.x + v.y*v.y;
      const mic_f u = clamp(d0*rcp(d1),zero,one);
      const mic4f p = p0 + u*v;
      const mic_f t = p.z * inv_ray_length;
      const mic_f d2 = p.x*p.x + p.y*p.y; 
      const mic_f r = p.w;
      const mic_f r2 = r*r;
      mic_m valid = le(d2,r2);
      valid = lt(valid,mic_f(ray.tnear),t);
      valid = lt(valid,t,mic_f(ray.tfar));


      if (unlikely(none(valid))) return false;
      STAT3(normal.trav_prim_hits,1,1,1);

      /* ray masking test */
#if defined(RTCORE_RAY_MASK)
      BezierCurves* g = ((Scene*)geom)->getBezierCurves(curve_in.geomID);
      if (unlikely(g->mask & ray.mask) == 0) return false;
#endif  

      if (ENABLE_INTERSECTION_FILTER) 
	{
	  const Geometry* const gg = ((Scene*)geom)->get(curve_in.geomID);
	  if (unlikely(gg->hasIntersectionFilter1())) 
	    {
	      while(any(valid)) 
		{
		  unsigned int i = select_min(valid,t);
		  float uu = (float(i)+u[i])*one_over_width; 
		  mic_f P,T;
		  eval(uu,p0123,P,T);
		  assert( T != mic_f::zero() );

		  if (runIntersectionFilter1(gg,ray,mic_f(uu),mic_f(0.0f),mic_f(t[i]),mic_f(T[0]),mic_f(T[1]),mic_f(T[2]),(mic_m)((unsigned int)1 << i),curve_in.geomID,curve_in.primID))
		    break;
		  valid ^= (unsigned int)1 << i;
		}
	      if (unlikely(none(valid))) return false;
	    }
	}

      unsigned int i = select_min(valid,t);
      float uu = (float(i)+u[i])*one_over_width; 
      mic_f P,T;
      eval(uu,p0123,P,T);
      assert( T != mic_f::zero() );

      ray.update((mic_m)1,mic_f(t[i]),mic_f(uu),mic_f::zero(),swAAAA(T),swBBBB(T),swCCCC(T),curve_in.geomID,curve_in.primID);

      return true;
    }

    static __forceinline bool occluded(const mic_f& pre_vx, 
				       const mic_f& pre_vy, 
				       const mic_f& pre_vz, 					
				       const mic_f& inv_ray_length,
				       const Ray& ray, 
				       const mic_f &dir_xyz,
				       const mic_f &org_xyz,
				       const Bezier1i& curve_in, 
				       const void* geom) 
    {
      STAT3(shadow.trav_prims,1,1,1);
      const mic_f zero = mic_f::zero();
      const mic_f one  = mic_f::one();

      prefetch<PFHINT_L1>(curve_in.p + 0);
      prefetch<PFHINT_L1>(curve_in.p + 3);

      const mic_f p0123 = uload16f((float*)curve_in.p);
      const mic_f p0123_org = p0123 - org_xyz;

      const mic_f p0123_2D = select(0x7777,pre_vx * swAAAA(p0123_org) + pre_vy * swBBBB(p0123_org) + pre_vz * swCCCC(p0123_org),p0123);


      const mic_f c0 = load16f(&coeff01[0]);
      const mic_f c1 = load16f(&coeff01[1]);
      const mic_f c2 = load16f(&coeff01[2]);
      const mic_f c3 = load16f(&coeff01[3]);

      const mic4f p0 = eval16(p0123_2D,c0,c1,c2,c3);

      const mic_f last_x = mic_f(p0123_2D[12 + 0]);
      const mic_f last_y = mic_f(p0123_2D[13 + 0]);
      const mic_f last_z = mic_f(p0123_2D[14 + 0]);
      const mic_f last_w = mic_f(p0123_2D[15 + 0]);

      const mic4f p1(align_shift_right<1>(last_x,p0[0]),  
       		     align_shift_right<1>(last_y,p0[1]), 
      		     align_shift_right<1>(last_z,p0[2]),  
       		     align_shift_right<1>(last_w,p0[3]));

      /* approximative intersection with cone */
      const mic4f v = p1-p0;
      const mic4f w = -p0;
      const mic_f d0 = w.x*v.x + w.y*v.y;
      const mic_f d1 = v.x*v.x + v.y*v.y;
      const mic_f u = clamp(d0*rcp(d1),zero,one);
      const mic4f p = p0 + u*v;
      const mic_f t = p.z * inv_ray_length;
      const mic_f d2 = p.x*p.x + p.y*p.y; 
      const mic_f r = p.w;
      const mic_f r2 = r*r;
      mic_m valid = le(d2,r2);
      valid = lt(valid,mic_f(ray.tnear),t);
      valid = lt(valid,t,mic_f(ray.tfar));


      if (unlikely(none(valid))) return false;

      STAT3(shadow.trav_prim_hits,1,1,1);

      /* ray masking test */
#if defined(RTCORE_RAY_MASK)
      BezierCurves* g = ((Scene*)geom)->getBezierCurves(curve_in.geomID);
      if (unlikely(g->mask & ray.mask) == 0) return false;
#endif  

      if (ENABLE_INTERSECTION_FILTER) 
	{
	  const Geometry* const gg = ((Scene*)geom)->get(curve_in.geomID);
	  if (unlikely(gg->hasOcclusionFilter1())) 
	    {
	      while(any(valid)) 
		{
		  unsigned int i = select_min(valid,t);
		  const float one_over_width = 1.0f/16.0f;
		  float uu = (float(i)+u[i])*one_over_width; 
		  mic_f P,T;
		  eval(uu,p0123,P,T);
		  assert( T != mic_f::zero() );

		  if (runOcclusionFilter1(gg,(Ray&)ray,mic_f(uu),mic_f(0.0f),mic_f(t[i]),mic_f(T[0]),mic_f(T[1]),mic_f(T[2]),(mic_m)((unsigned int)1 << i),curve_in.geomID,curve_in.primID))
		    return true;
		  valid ^= (unsigned int)1 << i;
		}
	    }
	  return false;
	}

      return true;
    }
  };
}
