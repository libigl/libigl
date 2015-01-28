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

#ifndef __EMBREE_ACCEL_TRIANGLE4V_INTERSECTOR1_PLUECKER_H__
#define __EMBREE_ACCEL_TRIANGLE4V_INTERSECTOR1_PLUECKER_H__

#include "triangle4v.h"
#include "../common/ray.h"

namespace embree
{
  struct Triangle4vIntersector1Pluecker
  {
    typedef Triangle4v Primitive;

    /*! Intersect a ray with the 4 triangles and updates the hit. */
    static __forceinline void intersect(Ray& ray, const Triangle4v& tri, void* geom)
    {
      /* calculate vertices relative to ray origin */
      STAT3(normal.trav_prims,1,1,1);
      const sse3f O = sse3f(ray.org);
      const sse3f D = sse3f(ray.dir);
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
      valid &= (T >= absDen*ssef(ray.tnear)) & (absDen*ssef(ray.tfar) >= T);
      if (unlikely(none(valid))) return;

        /* perform backface culling */
#if defined(__BACKFACE_CULLING__)
      valid &= den > ssef(zero);
      if (unlikely(none(valid))) return;
#else
      valid &= den != ssef(zero);
      if (unlikely(none(valid))) return;
#endif

      /* ray masking test */
#if defined(__USE_RAY_MASK__)
      valid &= (tri.mask & ray.mask) != 0;
      if (unlikely(none(valid))) return;
#endif

      /* calculate hit information */
      const ssef u = U / absDen;
      const ssef v = V / absDen;
      const ssef t = T / absDen;
      size_t i = select_min(valid,t);
      int geomID = tri.geomID[i];
      
      /* intersection filter test */
#if defined(__INTERSECTION_FILTER__)
      while (true) 
      {
        Geometry* geometry = ((Scene*)geom)->get(geomID);
        if (likely(!geometry->hasIntersectionFilter1())) 
        {
#endif
          /* update hit information */
          ray.u = u[i];
          ray.v = v[i];
          ray.tfar = t[i];
          ray.Ng.x = Ng.x[i];
          ray.Ng.y = Ng.y[i];
          ray.Ng.z = Ng.z[i];
          ray.geomID = geomID;
          ray.primID = tri.primID[i];

#if defined(__INTERSECTION_FILTER__)
          return;
        }

        Vec3fa N = Vec3fa(Ng.x[i],Ng.y[i],Ng.z[i]);
        if (runIntersectionFilter1(geometry,ray,u[i],v[i],t[i],N,geomID,tri.primID[i])) return;
        valid[i] = 0;
        if (none(valid)) return;
        i = select_min(valid,t);
        geomID = tri.geomID[i];
      }
#endif
    }

    static __forceinline void intersect(Ray& ray, const Triangle4v* tri, size_t num, void* geom)
    {
      for (size_t i=0; i<num; i++)
        intersect(ray,tri[i],geom);
    }

    /*! Test if the ray is occluded by one of the triangles. */
    static __forceinline bool occluded(Ray& ray, const Triangle4v& tri, void* geom)
    {
      /* calculate vertices relative to ray origin */
      STAT3(shadow.trav_prims,1,1,1);
      const sse3f O = sse3f(ray.org);
      const sse3f D = sse3f(ray.dir);
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
      valid &= (T >= absDen*ssef(ray.tnear)) & (absDen*ssef(ray.tfar) >= T);
      if (unlikely(none(valid))) return false;

      /* perform backface culling */
#if defined(__BACKFACE_CULLING__)
      valid &= den > ssef(zero);
      if (unlikely(none(valid))) return false;
#else
      valid &= den != ssef(zero);
      if (unlikely(none(valid))) return false;
#endif

      /* ray masking test */
#if defined(__USE_RAY_MASK__)
      valid &= (tri.mask & ray.mask) != 0;
      if (unlikely(none(valid))) return false;
#endif

      /* intersection filter test */
#if defined(__INTERSECTION_FILTER__)
      size_t m=movemask(valid), i=__bsf(m);
      while (true)
      {  
        const int geomID = tri.geomID[i];
        Geometry* geometry = ((Scene*)geom)->get(geomID);

        /* if we have no filter then the test passes */
        if (likely(!geometry->hasOcclusionFilter1()))
          break;

        /* calculate hit information */
        const ssef rcpAbsDen = rcp(absDen);
        const ssef u = U * rcpAbsDen;
        const ssef v = V * rcpAbsDen;
        const ssef t = T * rcpAbsDen;
        const Vec3fa N = Vec3fa(Ng.x[i],Ng.y[i],Ng.z[i]);
        if (runOcclusionFilter1(geometry,ray,u[i],v[i],t[i],N,geomID,tri.primID[i])) 
          break;

        /* test if one more triangle hit */
        m=__btc(m,i); i=__bsf(m);
        if (m == 0) return false;
      }
#endif

      return true;
    }

    static __forceinline bool occluded(Ray& ray, const Triangle4v* tri, size_t num, void* geom) 
    {
      for (size_t i=0; i<num; i++) 
        if (occluded(ray,tri[i],geom))
          return true;

      return false;
    }
  };
}

#endif


