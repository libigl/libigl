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

#include "triangle8.h"
#include "common/ray.h"
#include "geometry/filter.h"

namespace embree
{
  namespace isa
  {
    /*! Intersector for a single ray with 8 triangles. This intersector
     *  implements a modified version of the Moeller Trumbore
     *  intersector from the paper "Fast, Minimum Storage Ray-Triangle
     *  Intersection". In contrast to the paper we precalculate some
     *  factors and factor the calculations differently to allow
     *  precalculating the cross product e1 x e2. The resulting
     *  algorithm is similar to the fastest one of the paper "Optimizing
     *  Ray-Triangle Intersection via Automated Search". */
    template<bool list>
      struct Triangle8Intersector1MoellerTrumbore
      {
        typedef Triangle8 Primitive;
        
        struct Precalculations {
          __forceinline Precalculations (const Ray& ray) {}
        };
        
        /*! Intersect a ray with the 8 triangles and updates the hit. */
        static __forceinline void intersect(const Precalculations& pre, Ray& ray, const Primitive& tri, Scene* scene)
        {
          /* calculate denominator */
          STAT3(normal.trav_prims,1,1,1);
          const avx3f O = avx3f(ray.org);
          const avx3f D = avx3f(ray.dir);
          const avx3f C = avx3f(tri.v0) - O;
          const avx3f R = cross(D,C);
          const avxf den = dot(avx3f(tri.Ng),D);
          const avxf absDen = abs(den);
          const avxf sgnDen = signmsk(den);
          
          /* perform edge tests */
          const avxf U = dot(R,avx3f(tri.e2)) ^ sgnDen;
          const avxf V = dot(R,avx3f(tri.e1)) ^ sgnDen;
          
          /* perform backface culling */
#if defined(RTCORE_BACKFACE_CULLING)
          avxb valid = (den > avxf(zero)) & (U >= 0.0f) & (V >= 0.0f) & (U+V<=absDen);
#else
          avxb valid = (den != avxf(zero)) & (U >= 0.0f) & (V >= 0.0f) & (U+V<=absDen);
#endif
          if (likely(none(valid))) return;
          
          /* perform depth test */
          const avxf T = dot(avx3f(tri.Ng),C) ^ sgnDen;
          valid &= (T > absDen*avxf(ray.tnear)) & (T < absDen*avxf(ray.tfar));
          if (likely(none(valid))) return;
          
          /* ray masking test */
#if defined(RTCORE_RAY_MASK)
          valid &= (tri.mask & ray.mask) != 0;
          if (unlikely(none(valid))) return;
#endif
          
          /* calculate hit information */
          const avxf rcpAbsDen = rcp(absDen);
          const avxf u = U * rcpAbsDen;
          const avxf v = V * rcpAbsDen;
          const avxf t = T * rcpAbsDen;
          size_t i = select_min(valid,t);
          int geomID = tri.geomID<list>(i);
          
          /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER)
          while (true) 
          {
            Geometry* geometry = scene->get(geomID);
            if (likely(!geometry->hasIntersectionFilter1())) 
            {
#endif
              /* update hit information */
              ray.u = u[i];
              ray.v = v[i];
              ray.tfar = t[i];
              ray.Ng.x = tri.Ng.x[i];
              ray.Ng.y = tri.Ng.y[i];
              ray.Ng.z = tri.Ng.z[i];
              ray.geomID = geomID;
              ray.primID = tri.primID<list>(i);
              
#if defined(RTCORE_INTERSECTION_FILTER)
              return;
            }
            
            Vec3fa Ng = Vec3fa(tri.Ng.x[i],tri.Ng.y[i],tri.Ng.z[i]);
            if (runIntersectionFilter1(geometry,ray,u[i],v[i],t[i],Ng,geomID,tri.primID<list>(i))) return;
            valid[i] = 0;
            if (none(valid)) return;
            i = select_min(valid,t);
            geomID = tri.geomID<list>(i);
          }
#endif
        }
        
        /*! Test if the ray is occluded by one of the triangles. */
        static __forceinline bool occluded(const Precalculations& pre, Ray& ray, const Primitive& tri, Scene* scene)
        {
          /* calculate denominator */
          STAT3(shadow.trav_prims,1,1,1);
          const avx3f O = avx3f(ray.org);
          const avx3f D = avx3f(ray.dir);
          const avx3f C = avx3f(tri.v0) - O;
          const avx3f R = cross(D,C);
          const avxf den = dot(avx3f(tri.Ng),D);
          const avxf absDen = abs(den);
          const avxf sgnDen = signmsk(den);
          
          /* perform edge tests */
          const avxf U = dot(R,avx3f(tri.e2)) ^ sgnDen;
          const avxf V = dot(R,avx3f(tri.e1)) ^ sgnDen;
          const avxf W = absDen-U-V;
          avxb valid = (U >= 0.0f) & (V >= 0.0f) & (W >= 0.0f);
          if (unlikely(none(valid))) return false;
          
          /* perform depth test */
          const avxf T = dot(avx3f(tri.Ng),C) ^ sgnDen;
          valid &= (T >= absDen*avxf(ray.tnear)) & (absDen*avxf(ray.tfar) >= T);
          if (unlikely(none(valid))) return false;
          
          /* perform backface culling */
#if defined(RTCORE_BACKFACE_CULLING)
          valid &= den > avxf(zero);
          if (unlikely(none(valid))) return false;
#else
          valid &= den != avxf(zero);
          if (unlikely(none(valid))) return false;
#endif
          
          /* ray masking test */
#if defined(RTCORE_RAY_MASK)
          valid &= (tri.mask & ray.mask) != 0;
          if (unlikely(none(valid))) return false;
#endif
          
          /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER)
          
          size_t i = select_min(valid,T);
          int geomID = tri.geomID<list>(i);
          
          while (true) 
          {
            Geometry* geometry = scene->get(geomID);
            if (likely(!geometry->hasOcclusionFilter1())) break;
            
            /* calculate hit information */
            const avxf rcpAbsDen = rcp(absDen);
            const avxf u = U * rcpAbsDen;
            const avxf v = V * rcpAbsDen;
            const avxf t = T * rcpAbsDen;
            const Vec3fa Ng = Vec3fa(tri.Ng.x[i],tri.Ng.y[i],tri.Ng.z[i]);
            if (runOcclusionFilter1(geometry,ray,u[i],v[i],t[i],Ng,geomID,tri.primID<list>(i))) break;
            valid[i] = 0;
            if (none(valid)) return false;
            i = select_min(valid,T);
            geomID = tri.geomID<list>(i);
          }
#endif
          
          return true;
        }
      };
  }
}
