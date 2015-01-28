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

#ifndef __EMBREE_ACCEL_TRIANGLE1_INTERSECTOR1V_MOELLER_MB_H__
#define __EMBREE_ACCEL_TRIANGLE1_INTERSECTOR1V_MOELLER_MB_H__

#include "triangle1v.h"
#include "common/ray.h"
#include "geometry/filter.h"

namespace embree
{
  /*! Intersector for a single ray with individual precomputed
   *  triangles. This intersector implements a modified version of the
   *  Moeller Trumbore intersector from the paper "Fast, Minimum
   *  Storage Ray-Triangle Intersection". In contrast to the paper we
   *  precalculate some factors and factor the calculations
   *  differently to allow precalculating the cross product e1 x
   *  e2. The resulting algorithm is similar to the fastest one of the
   *  paper "Optimizing Ray-Triangle Intersection via Automated
   *  Search". */
  struct Triangle1vIntersector1MoellerTrumboreMB
  {
    typedef Triangle1vMB Primitive;

    /*! Intersect a ray with the triangle and updates the hit. */
    static __forceinline void intersect(Ray& ray, const Triangle1vMB& tri, const void* geom)
    {
      /* load triangle */
      STAT3(normal.trav_prims,1,1,1);
      const Vec3fa tri_v0 = tri.v0+ray.time*tri.d0;
      const Vec3fa tri_v1 = tri.v1+ray.time*tri.d1;
      const Vec3fa tri_v2 = tri.v2+ray.time*tri.d2;
      const Vec3fa e1 = tri_v0-tri_v1;
      const Vec3fa e2 = tri_v2-tri_v0;
      const Vec3fa tri_Ng = cross(e1,e2);

      /* calculate denominator */
      const Vec3fa O = ray.org;
      const Vec3fa D = ray.dir;
      const Vec3fa C = tri_v0 - O;
      const Vec3fa R = cross(D,C);
      const float den = dot(tri_Ng,D);
      const float absDen = abs(den);
      const float sgnDen = signmsk(den);

      /* perform edge tests */
      const float U = xorf(dot(R,e2),sgnDen);
      if (unlikely(U < 0.0f)) return;
      const float V = xorf(dot(R,e1),sgnDen);
      if (unlikely(V < 0.0f)) return;
      const float W = absDen-U-V;
      if (unlikely(W < 0.0f)) return;
      
      /* perform depth test */
      const float T = xorf(dot(tri_Ng,C),sgnDen);
      if (unlikely(absDen*ray.tfar < T)) return;
      if (unlikely(T < absDen*ray.tnear)) return;

      /* perform backface culling */
#if defined(__BACKFACE_CULLING__)
      if (unlikely(den <= 0.0f)) return;
#else
      if (unlikely(den == 0.0f)) return;
#endif

      /* ray masking test */
#if defined(__USE_RAY_MASK__)
      if (unlikely((tri.mask() & ray.mask) == 0)) return;
#endif

      /* calculate hit information */
      const float rcpAbsDen = rcp(absDen);
      const float u = U * rcpAbsDen;
      const float v = V * rcpAbsDen;
      const float t = T * rcpAbsDen;
      const int geomID = tri.geomID();
      const int primID = tri.primID();

      /* intersection filter test */
#if defined(__INTERSECTION_FILTER__)
      Geometry* geometry = ((Scene*)geom)->get(geomID);
      if (unlikely(geometry->hasIntersectionFilter1())) {
        runIntersectionFilter1(geometry,ray,u,v,t,tri_Ng,geomID,primID);
        return;
      }
#endif

      /* update hit information */
      ray.u = u;
      ray.v = v;
      ray.tfar = t;
      ray.Ng  = tri_Ng;
      ray.geomID = geomID;
      ray.primID = primID;
    }

    static __forceinline void intersect(Ray& ray, const Triangle1vMB* tri, size_t num, void* geom)
    {
      for (size_t i=0; i<num; i++)
        intersect(ray,tri[i],geom);
    }

    /*! Test if the ray is occluded by one of the triangles. */
    static __forceinline bool occluded(Ray& ray, const Triangle1vMB& tri, const void* geom)
    {
      /* load triangle */
      STAT3(shadow.trav_prims,1,1,1);
      const Vec3fa tri_v0 = tri.v0+ray.time*tri.d0;
      const Vec3fa tri_v1 = tri.v1+ray.time*tri.d1;
      const Vec3fa tri_v2 = tri.v2+ray.time*tri.d2;
      const Vec3fa e1 = tri_v0-tri_v1;
      const Vec3fa e2 = tri_v2-tri_v0;
      const Vec3fa tri_Ng = cross(e1,e2);

      /* calculate denominator */
      const Vec3fa O = Vec3fa(ray.org);
      const Vec3fa D = Vec3fa(ray.dir);
      const Vec3fa C = tri_v0 - O;
      const Vec3fa R = cross(D,C);
      const float den = dot(tri_Ng,D);
      const float absDen = abs(den);
      const float sgnDen = signmsk(den);

      /* perform edge tests */
      const float U = xorf(dot(R,e2),sgnDen);
      if (unlikely(U < 0.0f)) return false;
      const float V = xorf(dot(R,e1),sgnDen);
      if (unlikely(V < 0.0f)) return false;
      const float W = absDen-U-V;
      if (unlikely(W < 0.0f)) return false;
      
      /* perform depth test */
      const float T = xorf(dot(tri_Ng,C),sgnDen);
      if (unlikely(absDen*ray.tfar < T)) return false;
      if (unlikely(T < absDen*ray.tnear)) return false;

      /* perform backface culling */
#if defined(__BACKFACE_CULLING__)
      if (unlikely(den <= 0.0f)) return false;
#else
      if (unlikely(den == 0.0f)) return false;
#endif

      /* ray masking test */
#if defined(__USE_RAY_MASK__)
      if (unlikely((tri.mask() & ray.mask) == 0)) return false;
#endif

      /* intersection filter test */
#if defined(__INTERSECTION_FILTER__)
      const int geomID = tri.geomID();
      Geometry* geometry = ((Scene*)geom)->get(geomID);
      if (unlikely(geometry->hasOcclusionFilter1()))
      {
        /* calculate hit information */
        const float rcpAbsDen = rcp(absDen);
        const float u = U*rcpAbsDen;
        const float v = V*rcpAbsDen;
        const float t = T*rcpAbsDen;
        const int primID = tri.primID();
        return runOcclusionFilter1(geometry,ray,u,v,t,tri_Ng,geomID,primID);
      }
#endif

      return true;
    }

    static __forceinline bool occluded(Ray& ray, const Triangle1vMB* tri, size_t num, void* geom) 
    {
      for (size_t i=0; i<num; i++) 
        if (occluded(ray,tri[i],geom))
          return true;

      return false;
    }
  };
}

#endif


