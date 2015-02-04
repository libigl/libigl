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
#include "triangle8_intersector1_moeller.h"
#include "../common/ray8.h"

namespace embree
{
  namespace isa
  {
    /*! Intersector for 8 triangles with 8 rays. This intersector
     *  implements a modified version of the Moeller Trumbore
     *  intersector from the paper "Fast, Minimum Storage Ray-Triangle
     *  Intersection". In contrast to the paper we precalculate some
     *  factors and factor the calculations differently to allow
     *  precalculating the cross product e1 x e2. */
    template<bool list, bool enableIntersectionFilter>
      struct Triangle8Intersector8MoellerTrumbore
      {
        typedef Triangle8 Primitive;
        
        struct Precalculations {
          __forceinline Precalculations (const avxb& valid, const Ray8& ray) {}
        };
        
        /*! Intersects a 8 rays with 8 triangles. */
        static __forceinline void intersect(const avxb& valid_i, Precalculations& pre, Ray8& ray, const Primitive& tri, Scene* scene)
        {
          for (size_t i=0; i<8; i++)
          {
            if (!tri.valid(i)) break;
            STAT3(normal.trav_prims,1,popcnt(valid_i),8);
            
            /* load edges and geometry normal */
            avxb valid = valid_i;
            const avx3f p0 = broadcast8f(tri.v0,i);
            const avx3f e1 = broadcast8f(tri.e1,i);
            const avx3f e2 = broadcast8f(tri.e2,i);
            const avx3f Ng = broadcast8f(tri.Ng,i);
            
            /* calculate denominator */
            const avx3f C = p0 - ray.org;
            const avx3f R = cross(ray.dir,C);
            const avxf den = dot(Ng,ray.dir);
            const avxf absDen = abs(den);
            const avxf sgnDen = signmsk(den);
            
            /* test against edge p2 p0 */
            const avxf U = dot(R,e2) ^ sgnDen;
            valid &= U >= 0.0f;
            
            /* test against edge p0 p1 */
            const avxf V = dot(R,e1) ^ sgnDen;
            valid &= V >= 0.0f;
            
            /* test against edge p1 p2 */
            const avxf W = absDen-U-V;
            valid &= W >= 0.0f;
            if (likely(none(valid))) continue;
            
            /* perform depth test */
            const avxf T = dot(Ng,C) ^ sgnDen;
            valid &= (T >= absDen*ray.tnear) & (absDen*ray.tfar >= T);
            if (unlikely(none(valid))) continue;
            
            /* perform backface culling */
#if defined(RTCORE_BACKFACE_CULLING)
            valid &= den > avxf(zero);
            if (unlikely(none(valid))) continue;
#else
            valid &= den != avxf(zero);
            if (unlikely(none(valid))) continue;
#endif
            
            /* ray masking test */
#if defined(RTCORE_RAY_MASK)
            valid &= (tri.mask[i] & ray.mask) != 0;
            if (unlikely(none(valid))) continue;
#endif
            
            /* calculate hit information */
            const avxf rcpAbsDen = rcp(absDen);
            const avxf u = U*rcpAbsDen;
            const avxf v = V*rcpAbsDen;
            const avxf t = T*rcpAbsDen;
            const int geomID = tri.geomID<list>(i);
            const int primID = tri.primID<list>(i);
            
            /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER)
            if (enableIntersectionFilter) {
              Geometry* geometry = scene->get(geomID);
              if (unlikely(geometry->hasIntersectionFilter8())) {
                runIntersectionFilter8(valid,geometry,ray,u,v,t,Ng,geomID,primID);
                continue;
              }
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
        
        /*! Test for 4 rays if they are occluded by any of the 4 triangle. */
        static __forceinline avxb occluded(const avxb& valid_i, Precalculations& pre, Ray8& ray, const Primitive& tri, Scene* scene)
        {
          avxb valid0 = valid_i;
          
          for (size_t i=0; i<8; i++)
          {
            if (!tri.valid(i)) break;
            STAT3(shadow.trav_prims,1,popcnt(valid_i),8);
            
            /* load edges and geometry normal */
            avxb valid = valid0;
            const avx3f p0 = broadcast8f(tri.v0,i);
            const avx3f e1 = broadcast8f(tri.e1,i);
            const avx3f e2 = broadcast8f(tri.e2,i);
            const avx3f Ng = broadcast8f(tri.Ng,i);
            
            /* calculate denominator */
            const avx3f C = p0 - ray.org;
            const avx3f R = cross(ray.dir,C);
            const avxf den = dot(Ng,ray.dir);
            const avxf absDen = abs(den);
            const avxf sgnDen = signmsk(den);
            
            /* test against edge p2 p0 */
            const avxf U = dot(R,e2) ^ sgnDen;
            valid &= U >= 0.0f;
            
            /* test against edge p0 p1 */
            const avxf V = dot(R,e1) ^ sgnDen;
            valid &= V >= 0.0f;
            
            /* test against edge p1 p2 */
            const avxf W = absDen-U-V;
            valid &= W >= 0.0f;
            if (likely(none(valid))) continue;
            
            /* perform depth test */
            const avxf T = dot(Ng,C) ^ sgnDen;
            valid &= (T >= absDen*ray.tnear) & (absDen*ray.tfar >= T);
            if (unlikely(none(valid))) continue;
            
            /* perform backface culling */
#if defined(RTCORE_BACKFACE_CULLING)
            valid &= den > avxf(zero);
            if (unlikely(none(valid))) continue;
#else
            valid &= den != avxf(zero);
            if (unlikely(none(valid))) continue;
#endif
            
            /* ray masking test */
#if defined(RTCORE_RAY_MASK)
            valid &= (tri.mask[i] & ray.mask) != 0;
            if (unlikely(none(valid))) continue;
#endif
            
            /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER)
            if (enableIntersectionFilter) 
            {
              const int geomID = tri.geomID<list>(i);
              Geometry* geometry = scene->get(geomID);
              if (unlikely(geometry->hasOcclusionFilter8()))
              {
                /* calculate hit information */
                const avxf rcpAbsDen = rcp(absDen);
                const avxf u = U*rcpAbsDen;
                const avxf v = V*rcpAbsDen;
                const avxf t = T*rcpAbsDen;
                const int primID = tri.primID<list>(i);
                valid = runOcclusionFilter8(valid,geometry,ray,u,v,t,Ng,geomID,primID);
              }
            }
#endif
            
            /* update occlusion */
            valid0 &= !valid;
            if (none(valid0)) break;
          }
          return !valid0;
        }
        
        /*! Intersect a ray with the 4 triangles and updates the hit. */
        static __forceinline void intersect(Precalculations& pre, Ray8& ray, size_t k, const Primitive& tri, Scene* scene)
        {
          /* calculate denominator */
          STAT3(normal.trav_prims,1,1,1);
          const avx3f O = broadcast8f(ray.org,k);
          const avx3f D = broadcast8f(ray.dir,k);
          const avx3f C = avx3f(tri.v0) - O;
          const avx3f R = cross(D,C);
          const avxf den = dot(avx3f(tri.Ng),D);
          const avxf absDen = abs(den);
          const avxf sgnDen = signmsk(den);
          
          /* perform edge tests */
          const avxf U = dot(R,avx3f(tri.e2)) ^ sgnDen;
          const avxf V = dot(R,avx3f(tri.e1)) ^ sgnDen;
          avxb valid = (U >= 0.0f) & (V >= 0.0f) & (U+V<=absDen);
          if (likely(none(valid))) return;
          
          /* perform depth test */
          const avxf T = dot(avx3f(tri.Ng),C) ^ sgnDen;
          valid &= (T > absDen*avxf(ray.tnear[k])) & (T < absDen*avxf(ray.tfar[k]));
          if (likely(none(valid))) return;
          
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
          valid &= (tri.mask & ray.mask[k]) != 0;
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
            if (likely(!enableIntersectionFilter || !geometry->hasIntersectionFilter8())) 
            {
#endif
              /* update hit information */
              ray.u[k] = u[i];
              ray.v[k] = v[i];
              ray.tfar[k] = t[i];
              ray.Ng.x[k] = tri.Ng.x[i];
              ray.Ng.y[k] = tri.Ng.y[i];
              ray.Ng.z[k] = tri.Ng.z[i];
              ray.geomID[k] = geomID;
              ray.primID[k] = tri.primID<list>(i);
              
#if defined(RTCORE_INTERSECTION_FILTER)
              return;
            }
            
            const Vec3fa Ng(tri.Ng.x[i],tri.Ng.y[i],tri.Ng.z[i]);
            if (runIntersectionFilter8(geometry,ray,k,u[i],v[i],t[i],Ng,geomID,tri.primID<list>(i))) return;
            valid[i] = 0;
            if (unlikely(none(valid))) return;
            i = select_min(valid,t);
            geomID = tri.geomID<list>(i);
          }
#endif
        }
        
        /*! Test if the ray is occluded by one of the triangles. */
        static __forceinline bool occluded(Precalculations& pre, Ray8& ray, size_t k, const Primitive& tri, Scene* scene)
        {
          /* calculate denominator */
          STAT3(shadow.trav_prims,1,1,1);
          const avx3f O = broadcast8f(ray.org,k);
          const avx3f D = broadcast8f(ray.dir,k);
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
          valid &= (T >= absDen*avxf(ray.tnear[k])) & (absDen*avxf(ray.tfar[k]) >= T);
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
          valid &= (tri.mask & ray.mask[k]) != 0;
          if (unlikely(none(valid))) return false;
#endif
          
          /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER)
          
          size_t i = select_min(valid,T);
          int geomID = tri.geomID<list>(i);
          
          while (true) 
          {
            Geometry* geometry = scene->get(geomID);
            if (likely(!enableIntersectionFilter || !geometry->hasOcclusionFilter8())) break;
            
            /* calculate hit information */
            const avxf rcpAbsDen = rcp(absDen);
            const avxf u = U * rcpAbsDen;
            const avxf v = V * rcpAbsDen;
            const avxf t = T * rcpAbsDen;
            const Vec3fa Ng(tri.Ng.x[i],tri.Ng.y[i],tri.Ng.z[i]);
            if (runOcclusionFilter8(geometry,ray,k,u[i],v[i],t[i],Ng,geomID,tri.primID<list>(i))) break;
            valid[i] = 0;
            if (unlikely(none(valid))) return false;
            i = select_min(valid,T);
            geomID = tri.geomID<list>(i);
          }
#endif
          
          return true;
        }
      };
  }
}
