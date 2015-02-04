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

#include "triangle4v.h"
#include "../common/ray4.h"

namespace embree
{
  namespace isa
  {
    template<bool list>
      struct Triangle4vIntersector4Pluecker
      {
        typedef Triangle4v Primitive;
        
        struct Precalculations {
          __forceinline Precalculations (const sseb& valid, const Ray4& ray) {}
        };
        
        /*! Intersects a 4 rays with 4 triangles. */
        static __forceinline void intersect(const sseb& valid_i, Precalculations& pre, Ray4& ray, const Primitive& tri, Scene* scene)
        {
          for (size_t i=0; i<4; i++)
          {
            if (!tri.valid(i)) break;
            STAT3(normal.trav_prims,1,popcnt(valid_i),4);
            
            /* calculate vertices relative to ray origin */
            sseb valid = valid_i;
            const sse3f O = ray.org;
            const sse3f D = ray.dir;
            const sse3f v0 = broadcast4f(tri.v0,i)-O;
            const sse3f v1 = broadcast4f(tri.v1,i)-O;
            const sse3f v2 = broadcast4f(tri.v2,i)-O;
            
            /* calculate triangle edges */
            const sse3f e0 = v2-v0;
            const sse3f e1 = v0-v1;
            const sse3f e2 = v1-v2;
            
            /* calculate geometry normal and denominator */
            const sse3f Ng1 = cross(e1,e0);
            const sse3f Ng = Ng1+Ng1;
            const ssef den = dot(sse3f(Ng),D);
            const ssef absDen = abs(den);
            const ssef sgnDen = signmsk(den);
            
            /* perform edge tests */
            const ssef U = dot(sse3f(cross(v2+v0,e0)),D) ^ sgnDen;
            valid &= U >= 0.0f;
            if (likely(none(valid))) continue;
            const ssef V = dot(sse3f(cross(v0+v1,e1)),D) ^ sgnDen;
            valid &= V >= 0.0f;
            if (likely(none(valid))) continue;
            const ssef W = dot(sse3f(cross(v1+v2,e2)),D) ^ sgnDen;
            valid &= W >= 0.0f;
            if (likely(none(valid))) continue;
            
            /* perform depth test */
            const ssef T = dot(v0,sse3f(Ng)) ^ sgnDen;
            valid &= (T >= absDen*ray.tnear) & (absDen*ray.tfar >= T);
            if (unlikely(none(valid))) continue;
            
            /* perform backface culling */
#if defined(RTCORE_BACKFACE_CULLING)
            valid &= den > ssef(zero);
            if (unlikely(none(valid))) continue;
#else
            valid &= den != ssef(zero);
            if (unlikely(none(valid))) continue;
#endif
            
            /* ray masking test */
#if defined(RTCORE_RAY_MASK)
            valid &= (tri.mask[i] & ray.mask) != 0;
            if (unlikely(none(valid))) continue;
#endif
            
            /* calculate hit information */
            const ssef u = U / absDen;
            const ssef v = V / absDen;
            const ssef t = T / absDen;
            const int geomID = tri.geomID<list>(i);
            const int primID = tri.primID<list>(i);
            
            /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER)
            Geometry* geometry = scene->get(geomID);
            if (unlikely(geometry->hasIntersectionFilter4())) {
              runIntersectionFilter4(valid,geometry,ray,u,v,t,Ng,geomID,primID);
              continue;
            }
#endif
            
            /* update hit information */
            store4f(valid,&ray.u,u);
            store4f(valid,&ray.v,v);
            store4f(valid,&ray.tfar,t);
            store4i(valid,&ray.geomID,geomID);
            store4i(valid,&ray.primID,primID);
            store4f(valid,&ray.Ng.x,Ng.x);
            store4f(valid,&ray.Ng.y,Ng.y);
            store4f(valid,&ray.Ng.z,Ng.z);
          }
        }
        
        /*! Test for 4 rays if they are occluded by any of the 4 triangle. */
        static __forceinline sseb occluded(const sseb& valid_i, Precalculations& pre, Ray4& ray, const Primitive& tri, Scene* scene)
        {
          sseb valid0 = valid_i;
          
          for (size_t i=0; i<4; i++)
          {
            if (!tri.valid(i)) break;
            STAT3(shadow.trav_prims,1,popcnt(valid_i),4);
            
            /* calculate vertices relative to ray origin */
            sseb valid = valid0;
            const sse3f O = ray.org;
            const sse3f D = ray.dir;
            const sse3f v0 = broadcast4f(tri.v0,i)-O;
            const sse3f v1 = broadcast4f(tri.v1,i)-O;
            const sse3f v2 = broadcast4f(tri.v2,i)-O;
            
            /* calculate triangle edges */
            const sse3f e0 = v2-v0;
            const sse3f e1 = v0-v1;
            const sse3f e2 = v1-v2;
            
            /* calculate geometry normal and denominator */
            const sse3f Ng1 = cross(e1,e0);
            const sse3f Ng = Ng1+Ng1;
            const ssef den = dot(sse3f(Ng),D);
            const ssef absDen = abs(den);
            const ssef sgnDen = signmsk(den);
            
            /* perform edge tests */
            const ssef U = dot(sse3f(cross(v2+v0,e0)),D) ^ sgnDen;
            valid &= U >= 0.0f;
            if (likely(none(valid))) continue;
            const ssef V = dot(sse3f(cross(v0+v1,e1)),D) ^ sgnDen;
            valid &= V >= 0.0f;
            if (likely(none(valid))) continue;
            const ssef W = dot(sse3f(cross(v1+v2,e2)),D) ^ sgnDen;
            valid &= W >= 0.0f;
            if (likely(none(valid))) continue;
            
            /* perform depth test */
            const ssef T = dot(v0,sse3f(Ng)) ^ sgnDen;
            valid &= (T >= absDen*ray.tnear) & (absDen*ray.tfar >= T);
            
            /* perform backface culling */
#if defined(RTCORE_BACKFACE_CULLING)
            valid &= den > ssef(zero);
            if (unlikely(none(valid))) continue;
#else
            valid &= den != ssef(zero);
            if (unlikely(none(valid))) continue;
#endif
            
            /* ray masking test */
#if defined(RTCORE_RAY_MASK)
            valid &= (tri.mask[i] & ray.mask) != 0;
            if (unlikely(none(valid))) continue;
#endif
            
            /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER)
            const int geomID = tri.geomID<list>(i);
            Geometry* geometry = scene->get(geomID);
            if (unlikely(geometry->hasOcclusionFilter4()))
            {
              /* calculate hit information */
              const ssef rcpAbsDen = rcp(absDen);
              const ssef u = U / absDen;
              const ssef v = V / absDen;
              const ssef t = T / absDen;
              const int primID = tri.primID<list>(i);
              valid = runOcclusionFilter4(valid,geometry,ray,u,v,t,Ng,geomID,primID);
            }
#endif
            
            /* update occlusion */
            valid0 &= !valid;
            if (none(valid0)) break;
          }
          return !valid0;
        }
        
        /*! Intersect a ray with the 4 triangles and updates the hit. */
        static __forceinline void intersect(Precalculations& pre, Ray4& ray, size_t k, const Primitive& tri, Scene* scene)
        {
          /* calculate vertices relative to ray origin */
          STAT3(normal.trav_prims,1,1,1);
          const sse3f O = broadcast4f(ray.org,k);
          const sse3f D = broadcast4f(ray.dir,k);
          const sse3f v0 = tri.v0-O;
          const sse3f v1 = tri.v1-O;
          const sse3f v2 = tri.v2-O;
          
          /* calculate triangle edges */
          const sse3f e0 = v2-v0;
          const sse3f e1 = v0-v1;
          const sse3f e2 = v1-v2;
          
          /* calculate geometry normal and denominator */
          const sse3f Ng1 = cross(e1,e0);
          const sse3f Ng = Ng1+Ng1;
          const ssef den = dot(Ng,D);
          const ssef absDen = abs(den);
          const ssef sgnDen = signmsk(den);
          
          /* perform edge tests */
          const ssef U = dot(cross(v2+v0,e0),D) ^ sgnDen;
          const ssef V = dot(cross(v0+v1,e1),D) ^ sgnDen;
          const ssef W = dot(cross(v1+v2,e2),D) ^ sgnDen;
          sseb valid = (U >= 0.0f) & (V >= 0.0f) & (W >= 0.0f);
          if (unlikely(none(valid))) return;
          
          /* perform depth test */
          const ssef T = dot(v0,Ng) ^ sgnDen;
          valid &= (T >= absDen*ssef(ray.tnear[k])) & (absDen*ssef(ray.tfar[k]) >= T);
          if (unlikely(none(valid))) return;
          
          /* perform backface culling */
#if defined(RTCORE_BACKFACE_CULLING)
          valid &= den > ssef(zero);
          if (unlikely(none(valid))) return;
#else
          valid &= den != ssef(zero);
          if (unlikely(none(valid))) return;
#endif
          
          /* ray masking test */
#if defined(RTCORE_RAY_MASK)
          valid &= (tri.mask & ray.mask[k]) != 0;
          if (unlikely(none(valid))) return;
#endif
          
          /* calculate hit information */
          const ssef u = U / absDen;
          const ssef v = V / absDen;
          const ssef t = T / absDen;
          size_t i = select_min(valid,t);
          int geomID = tri.geomID<list>(i);
          
          /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER)
          while (true) 
          {
            Geometry* geometry = scene->get(geomID);
            if (likely(!geometry->hasIntersectionFilter4())) 
            {
#endif
              /* update hit information */
              ray.u[k] = u[i];
              ray.v[k] = v[i];
              ray.tfar[k] = t[i];
              ray.Ng.x[k] = Ng.x[i];
              ray.Ng.y[k] = Ng.y[i];
              ray.Ng.z[k] = Ng.z[i];
              ray.geomID[k] = geomID;
              ray.primID[k] = tri.primID<list>(i);
              
#if defined(RTCORE_INTERSECTION_FILTER)
              return;
            }
            
            const Vec3fa N(Ng.x[i],Ng.y[i],Ng.z[i]);
            if (runIntersectionFilter4(geometry,ray,k,u[i],v[i],t[i],N,geomID,tri.primID<list>(i))) return;
            valid[i] = 0;
            if (unlikely(none(valid))) return;
            i = select_min(valid,t);
            geomID = tri.geomID<list>(i);
          }
#endif
        }
        
        /*! Test if the ray is occluded by one of the triangles. */
        static __forceinline bool occluded(Precalculations& pre, Ray4& ray, size_t k, const Primitive& tri, Scene* scene)
        {
          /* calculate vertices relative to ray origin */
          STAT3(shadow.trav_prims,1,1,1);
          const sse3f O = broadcast4f(ray.org,k);
          const sse3f D = broadcast4f(ray.dir,k);
          const sse3f v0 = tri.v0-O;
          const sse3f v1 = tri.v1-O;
          const sse3f v2 = tri.v2-O;
          
          /* calculate triangle edges */
          const sse3f e0 = v2-v0;
          const sse3f e1 = v0-v1;
          const sse3f e2 = v1-v2;
          
          /* calculate geometry normal and denominator */
          const sse3f Ng1 = cross(e1,e0);
          const sse3f Ng = Ng1+Ng1;
          const ssef den = dot(Ng,D);
          const ssef absDen = abs(den);
          const ssef sgnDen = signmsk(den);
          
          /* perform edge tests */
          const ssef U = dot(cross(v2+v0,e0),D) ^ sgnDen;
          const ssef V = dot(cross(v0+v1,e1),D) ^ sgnDen;
          const ssef W = dot(cross(v1+v2,e2),D) ^ sgnDen;
          sseb valid = (U >= 0.0f) & (V >= 0.0f) & (W >= 0.0f);
          if (unlikely(none(valid))) return false;
          
          /* perform depth test */
          const ssef T = dot(v0,Ng) ^ sgnDen;
          valid &= (T >= absDen*ssef(ray.tnear[k])) & (absDen*ssef(ray.tfar[k]) >= T);
          if (unlikely(none(valid))) return false;
          
          /* perform backface culling */
#if defined(RTCORE_BACKFACE_CULLING)
          valid &= den > ssef(zero);
          if (unlikely(none(valid))) return false;
#else
          valid &= den != ssef(zero);
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
            if (likely(!geometry->hasOcclusionFilter4())) break;
            
            /* calculate hit information */
            const ssef u = U / absDen;
            const ssef v = V / absDen;
            const ssef t = T / absDen;
            const Vec3fa N(Ng.x[i],Ng.y[i],Ng.z[i]);
            if (runOcclusionFilter4(geometry,ray,k,u[i],v[i],t[i],N,geomID,tri.primID<list>(i))) break;
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
