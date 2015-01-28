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

#ifndef __EMBREE_ACCEL_TRIANGLE1V_INTERSECTOR8_MOELLER_MB_H__
#define __EMBREE_ACCEL_TRIANGLE1V_INTERSECTOR8_MOELLER_MB_H__

#include "triangle1v.h"
#include "common/ray8.h"
#include "geometry/filter.h"

namespace embree
{
  /*! Intersector for individual precomputed triangles with 8
   *  rays. This intersector implements a modified version of the
   *  Moeller Trumbore intersector from the paper "Fast, Minimum
   *  Storage Ray-Triangle Intersection". In contrast to the paper we
   *  precalculate some factors and factor the calculations
   *  differently to allow precalculating the cross product e1 x
   *  e2. */
  struct Triangle1vIntersector8MoellerTrumboreMB
  {
    typedef Triangle1vMB Primitive;
    
    static __forceinline void intersect(const avxb& valid_i, Ray8& ray, const Triangle1vMB* __restrict__ tris, size_t num, const void* geom)
    {
      for (size_t i=0; i<num; i++) 
      {
        STAT3(normal.trav_prims,1,popcnt(valid_i),8);

        avxb valid = valid_i;
        const avx3f org = ray.org;
        const avx3f dir = ray.dir;
        const avxf zero = 0.0f;
        const Triangle1vMB& tri = tris[i];
        
        /* load vertices and calculate edges */
        const avx3f v0 = avx3f(tri.v0)+ray.time*avx3f(tri.d0);
        const avx3f v1 = avx3f(tri.v1)+ray.time*avx3f(tri.d1);
        const avx3f v2 = avx3f(tri.v2)+ray.time*avx3f(tri.d2);
        const avx3f e1 = v0-v1;
        const avx3f e2 = v2-v0;
        
        /* calculate denominator */
        const avx3f C =  v0 - org;
        const avx3f Ng = cross(e1,e2);
        const avxf den = dot(dir,Ng);
        const avxf sgnDen = signmsk(den);
        const avxf absDen = abs(den);
        
        /* perform edge tests */
        const avx3f R = cross(dir,C);
        const avxf V = dot(R,e1)^sgnDen;
        const avxf U = dot(R,e2)^sgnDen;
        valid &= V >= zero & U >= zero & U+V <= absDen;
        if (unlikely(none(valid))) continue;
      
        /* perform depth test */
        const avxf T = dot(C,Ng) ^ sgnDen;
        valid &= (T >= absDen*ray.tnear) & (absDen*ray.tfar >= T);
        if (unlikely(none(valid))) continue;

        /* perform backface culling */
#if defined(__BACKFACE_CULLING__)
        valid &= den > avxf(zero);
        if (unlikely(none(valid))) continue;
#else
        valid &= den != avxf(zero);
        if (unlikely(none(valid))) continue;
#endif
        
        /* ray masking test */
#if defined(__USE_RAY_MASK__)
        valid &= (tri.mask() & ray.mask) != 0;
        if (unlikely(none(valid))) continue;
#endif

        /* calculate hit information */
        const avxf rcpAbsDen = rcp(absDen);
        const avxf u = U*rcpAbsDen;
        const avxf v = V*rcpAbsDen;
        const avxf t = T*rcpAbsDen;
        const int geomID = tri.geomID();
        const int primID = tri.primID();

        /* intersection filter test */
#if defined(__INTERSECTION_FILTER__)
        Geometry* geometry = ((Scene*)geom)->get(geomID);
        if (unlikely(geometry->hasIntersectionFilter8())) {
          runIntersectionFilter8(valid,geometry,ray,u,v,t,Ng,geomID,primID);
          continue;
        }
#endif

        /* update hit information */
        store8f(valid,&ray.u,u);
        store8f(valid,&ray.v,v);
        store8f(valid,&ray.tfar,t);
        store8i(valid,&ray.geomID,geomID);
        store8i(valid,&ray.primID,primID);
        store8f(valid,&ray.Ng.x,Ng.x);
        store8f(valid,&ray.Ng.y,Ng.y);
        store8f(valid,&ray.Ng.z,Ng.z);
      }
    }

    static __forceinline avxb occluded(const avxb& valid_i, Ray8& ray, const Triangle1vMB* __restrict__ tris, size_t num, const void* geom)
    {
      avxb valid0 = valid_i;

      for (size_t i=0; i<num; i++) 
      {
        STAT3(shadow.trav_prims,1,popcnt(valid0),8);

        avxb valid = valid0;
        const avx3f org = ray.org;
        const avx3f dir = ray.dir;
        const avxf zero = 0.0f;
        const Triangle1vMB& tri = tris[i];
        
        /* load vertices and calculate edges */
        const avx3f v0 = avx3f(tri.v0)+ray.time*avx3f(tri.d0);
        const avx3f v1 = avx3f(tri.v1)+ray.time*avx3f(tri.d1);
        const avx3f v2 = avx3f(tri.v2)+ray.time*avx3f(tri.d2);
        const avx3f e1 = v0-v1;
        const avx3f e2 = v2-v0;
        
        /* calculate denominator */
        const avx3f C =  v0 - org;
        const avx3f Ng = cross(e1,e2);
        const avxf den = dot(dir,Ng);
        const avxf sgnDen = signmsk(den);
        const avxf absDen = abs(den);
        
        /* perform edge tests */
        const avx3f R = cross(dir,C);
        const avxf V = dot(R,e1)^sgnDen;
        const avxf U = dot(R,e2)^sgnDen;
        valid &= V >= zero & U >= zero & U+V <= absDen;
        if (unlikely(none(valid))) continue;
      
        /* perform depth test */
        const avxf T = dot(C,Ng) ^ sgnDen;
        valid &= (T >= absDen*ray.tnear) & (absDen*ray.tfar >= T);
        if (unlikely(none(valid))) continue;

        /* perform backface culling */
#if defined(__BACKFACE_CULLING__)
        valid &= den > avxf(zero);
        if (unlikely(none(valid))) continue;
#else
        valid &= den != avxf(zero);
        if (unlikely(none(valid))) continue;
#endif

        /* ray masking test */
#if defined(__USE_RAY_MASK__)
        valid &= (tri.mask() & ray.mask) != 0;
        if (unlikely(none(valid))) continue;
#endif

        /* intersection filter test */
#if defined(__INTERSECTION_FILTER__)
        const int geomID = tri.geomID();
        Geometry* geometry = ((Scene*)geom)->get(geomID);
        if (unlikely(geometry->hasOcclusionFilter8()))
        {
          /* calculate hit information */
          const avxf rcpAbsDen = rcp(absDen);
          const avxf u = U*rcpAbsDen;
          const avxf v = V*rcpAbsDen;
          const avxf t = T*rcpAbsDen;
          const int primID = tri.primID();
          valid = runOcclusionFilter8(valid,geometry,ray,u,v,t,Ng,geomID,primID);
        }
#endif
        
        /* update occlusion */
        valid0 &= !valid;
        if (none(valid0)) break;
      }
      return !valid0;
    }
  };
}

#endif


