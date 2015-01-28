// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
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

#ifndef __EMBREE_ACCEL_TRIANGLE1_INTERSECTOR16_MOELLER_H__
#define __EMBREE_ACCEL_TRIANGLE1_INTERSECTOR16_MOELLER_H__

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
  struct Triangle1Intersector16MoellerTrumbore
  {
    typedef Triangle1 Primitive;

    static __forceinline void intersect(const mic_m& valid_i, Ray16& ray, const Triangle1* __restrict__ tris, size_t num, const void* geom)
    {
      prefetch<PFHINT_L1>((mic_f*)tris +  0); 
      prefetch<PFHINT_L2>((mic_f*)tris +  1); 
      prefetch<PFHINT_L2>((mic_f*)tris +  2); 
      prefetch<PFHINT_L2>((mic_f*)tris +  3); 

      const mic_f zero = mic_f::zero();
      const mic_f one  = mic_f::one();

      const mic3f org = ray.org;
      const mic3f dir = ray.dir;

      for (size_t i=0; i<num; i++) 
      {
        const Triangle1& tri = tris[i];

        prefetch<PFHINT_L1>(&tris[i+1]); 

        STAT3(normal.trav_prims,1,popcnt(valid_i),16);

        mic_m valid = valid_i;
        
        /* load vertices and calculate edges */
        const mic_f v0 = broadcast4to16f(&tri.v0);
        const mic_f v1 = broadcast4to16f(&tri.v1);
        const mic_f v2 = broadcast4to16f(&tri.v2);
        const mic_f e1 = v0-v1;
        const mic_f e2 = v2-v0;

        /* calculate denominator */
        const mic3f _v0 = mic3f(swizzle<0>(v0),swizzle<1>(v0),swizzle<2>(v0));
        const mic3f C =  _v0 - org;
        const mic3f Ng = mic3f(tri.Ng);
        const mic_f den = dot(dir,Ng);

#if defined(__BACKFACE_CULLING__)
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
	prefetch<PFHINT_L1EX>(&ray.u);      
	prefetch<PFHINT_L1EX>(&ray.v);      
	prefetch<PFHINT_L1EX>(&ray.tfar);      

        if (unlikely(none(valid))) continue;
        const mic_f t = dot(C,Ng) * rcp_den;
      
        /* perform depth test */
        valid = ge(valid, t,ray.tnear);
	valid = ge(valid,ray.tfar,t);

	const mic_i geomID = tri.geomID();
	const mic_i primID = tri.primID();
	prefetch<PFHINT_L1EX>(&ray.geomID);      
	prefetch<PFHINT_L1EX>(&ray.primID);      
	prefetch<PFHINT_L1EX>(&ray.Ng.x);      
	prefetch<PFHINT_L1EX>(&ray.Ng.y);      
	prefetch<PFHINT_L1EX>(&ray.Ng.z);      

        /* ray masking test */
#if defined(__USE_RAY_MASK__)
        valid &= (tri.mask() & ray.mask) != 0;
#endif
        if (unlikely(none(valid))) continue;

        /* intersection filter test */
#if defined(__INTERSECTION_FILTER__)
        Geometry* geometry = ((Scene*)geom)->get(tri.geomID());
        if (unlikely(geometry->hasIntersectionFilter16())) {
          runIntersectionFilter16(valid,geometry,ray,u,v,t,Ng,geomID,primID);
          continue;
        }
#endif

        /* update hit information */
        store16f(valid,&ray.u,u);
        store16f(valid,&ray.v,v);
        store16f(valid,&ray.tfar,t);
        store16i(valid,&ray.geomID,geomID);
        store16i(valid,&ray.primID,primID);
        store16f(valid,&ray.Ng.x,Ng.x);
        store16f(valid,&ray.Ng.y,Ng.y);
        store16f(valid,&ray.Ng.z,Ng.z);
      }
    }

    static __forceinline mic_m occluded(const mic_m& valid_i, Ray16& ray, const Triangle1* __restrict__ tris, size_t num, const void* geom)
    {
      //prefetch<PFHINT_L1>((mic_f*)tris +  0); 
      prefetch<PFHINT_L2>((mic_f*)tris +  1); 
      prefetch<PFHINT_L2>((mic_f*)tris +  2); 
      prefetch<PFHINT_L2>((mic_f*)tris +  3); 

      mic_m valid0 = valid_i;
      const mic_f zero = mic_f::zero();
      const mic_f one  = mic_f::one();

      const mic3f org = ray.org;
      const mic3f dir = ray.dir;
     
      for (size_t i=0; i<num; i++) 
      {
        STAT3(shadow.trav_prims,1,popcnt(valid0),16);

        mic_m valid = valid0;
        const Triangle1& tri = tris[i];
        
        /* load vertices and calculate edges */
        const mic_f v0 = broadcast4to16f(&tri.v0);
        const mic_f v1 = broadcast4to16f(&tri.v1);
        const mic_f v2 = broadcast4to16f(&tri.v2);
        const mic_f e1 = v0-v1;
        const mic_f e2 = v2-v0;
        
        /* calculate denominator */
        const mic3f _v0 = mic3f(swizzle<0>(v0),swizzle<1>(v0),swizzle<2>(v0));
        const mic3f C =  _v0 - org;
        const mic_f Ng = broadcast4to16f(&tri.Ng);
        const mic3f _Ng = mic3f(swizzle<0>(Ng),swizzle<1>(Ng),swizzle<2>(Ng));
        const mic_f den = dot(dir,_Ng);

#if defined(__BACKFACE_CULLING__)
        valid &= den > zero;
#endif
        const mic_f rcp_den = rcp(den);
        const mic3f R = cross(dir,C);
        const mic3f _e1(swizzle<0>(e1),swizzle<1>(e1),swizzle<2>(e1));
        const mic_f u = dot(R,_e1)*rcp_den;
        const mic3f _e2(swizzle<0>(e2),swizzle<1>(e2),swizzle<2>(e2));
        const mic_f v = dot(R,_e2)*rcp_den;
	valid = ge(valid,u,zero);
	valid = ge(valid,v,zero);
	valid = le(valid,u+v,one); 
        const mic_f t = dot(C,_Ng) * rcp_den;

        if (unlikely(none(valid))) continue;
      
        /* perform depth test */
        valid = ge(valid, t,ray.tnear);
	valid = ge(valid,ray.tfar,t);

        /* ray masking test */
#if defined(__USE_RAY_MASK__)
        valid &= (tri.mask() & ray.mask) != 0;
#endif
        if (unlikely(none(valid))) continue;

        /* intersection filter test */
#if defined(__INTERSECTION_FILTER__)
        const mic_i geomID = tri.geomID();
	const mic_i primID = tri.primID();
        Geometry* geometry = ((Scene*)geom)->get(tri.geomID());
        if (unlikely(geometry->hasOcclusionFilter16()))
          valid = runOcclusionFilter16(valid,geometry,ray,u,v,t,Ng,geomID,primID);
#endif

        /* update occlusion */
        valid0 &= !valid;
        if (unlikely(none(valid0))) break;
      }
      return !valid0;
    }
  };
}

#endif


