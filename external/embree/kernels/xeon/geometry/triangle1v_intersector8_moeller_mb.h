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

#include "triangle1v.h"
#include "common/ray8.h"
#include "geometry/filter.h"

namespace embree
{
  namespace isa
  {
    /*! Intersector for individual precomputed triangles with 8
     *  rays. This intersector implements a modified version of the
     *  Moeller Trumbore intersector from the paper "Fast, Minimum
     *  Storage Ray-Triangle Intersection". In contrast to the paper we
     *  precalculate some factors and factor the calculations
     *  differently to allow precalculating the cross product e1 x
     *  e2. */
    template<bool list>
      struct Triangle1vIntersector8MoellerTrumboreMB
      {
        typedef Triangle1vMB Primitive;
        
        struct Precalculations {
          __forceinline Precalculations (const avxb& valid, const Ray8& ray) {}
        };
        
        static __forceinline void intersect(const avxb& valid_i, Precalculations& pre, Ray8& ray, const Primitive& tri, Scene* scene)
        {
          STAT3(normal.trav_prims,1,popcnt(valid_i),8);
          
          avxb valid = valid_i;
          const avx3f org = ray.org;
          const avx3f dir = ray.dir;
          const avxf zero = 0.0f;
          
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
          if (unlikely(none(valid))) return;
          
          /* perform depth test */
          const avxf T = dot(C,Ng) ^ sgnDen;
          valid &= (T >= absDen*ray.tnear) & (absDen*ray.tfar >= T);
          if (unlikely(none(valid))) return;
          
          /* perform backface culling */
#if defined(RTCORE_BACKFACE_CULLING)
          valid &= den > avxf(zero);
          if (unlikely(none(valid))) return;
#else
          valid &= den != avxf(zero);
          if (unlikely(none(valid))) return;
#endif
          
          /* ray masking test */
#if defined(RTCORE_RAY_MASK)
          valid &= (tri.mask() & ray.mask) != 0;
          if (unlikely(none(valid))) return;
#endif
          
          /* calculate hit information */
          const avxf rcpAbsDen = rcp(absDen);
          const avxf u = U*rcpAbsDen;
          const avxf v = V*rcpAbsDen;
          const avxf t = T*rcpAbsDen;
          const int geomID = tri.geomID<list>();
          const int primID = tri.primID<list>();
          
          /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER)
          Geometry* geometry = scene->get(geomID);
          if (unlikely(geometry->hasIntersectionFilter8())) {
            runIntersectionFilter8(valid,geometry,ray,u,v,t,Ng,geomID,primID);
            return;
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
        
        static __forceinline avxb occluded(const avxb& valid_i, Precalculations& pre, Ray8& ray, const Primitive& tri, Scene* scene)
        {
          STAT3(shadow.trav_prims,1,popcnt(valid_i),8);
          
          avxb valid = valid_i;
          const avx3f org = ray.org;
          const avx3f dir = ray.dir;
          const avxf zero = 0.0f;
          
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
          if (unlikely(none(valid))) return valid;
          
          /* perform depth test */
          const avxf T = dot(C,Ng) ^ sgnDen;
          valid &= (T >= absDen*ray.tnear) & (absDen*ray.tfar >= T);
          if (unlikely(none(valid))) return valid;
          
          /* perform backface culling */
#if defined(RTCORE_BACKFACE_CULLING)
          valid &= den > avxf(zero);
          if (unlikely(none(valid))) return valid;
#else
          valid &= den != avxf(zero);
          if (unlikely(none(valid))) return valid;
#endif
          
          /* ray masking test */
#if defined(RTCORE_RAY_MASK)
          valid &= (tri.mask() & ray.mask) != 0;
          if (unlikely(none(valid))) return valid;
#endif
          
          /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER)
          const int geomID = tri.geomID<list>();
          Geometry* geometry = scene->get(geomID);
          if (unlikely(geometry->hasOcclusionFilter8()))
          {
            /* calculate hit information */
            const avxf rcpAbsDen = rcp(absDen);
            const avxf u = U*rcpAbsDen;
            const avxf v = V*rcpAbsDen;
            const avxf t = T*rcpAbsDen;
            const int primID = tri.primID<list>();
            valid = runOcclusionFilter8(valid,geometry,ray,u,v,t,Ng,geomID,primID);
          }
#endif
          return valid;
        }
      };
  }
}
