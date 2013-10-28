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

#ifndef __EMBREE_ACCEL_TRIANGLE4_INTERSECTOR1_MOELLER_H__
#define __EMBREE_ACCEL_TRIANGLE4_INTERSECTOR1_MOELLER_H__

#include "triangle4.h"
#include "../common/ray.h"

namespace embree
{
  /*! Intersector for a single ray with 4 triangles. This intersector
   *  implements a modified version of the Moeller Trumbore
   *  intersector from the paper "Fast, Minimum Storage Ray-Triangle
   *  Intersection". In contrast to the paper we precalculate some
   *  factors and factor the calculations differently to allow
   *  precalculating the cross product e1 x e2. The resulting
   *  algorithm is similar to the fastest one of the paper "Optimizing
   *  Ray-Triangle Intersection via Automated Search". */
  struct Triangle4Intersector1MoellerTrumbore
  {
    typedef Triangle4 Triangle;

    /*! Name of intersector */
    static const char* name() { return "moeller"; }

    /*! Intersect a ray with the 4 triangles and updates the hit. */
    static __forceinline void intersect(Ray& ray, const Triangle4& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,1,1);

      /* calculate determinant */
      const sse3f O = sse3f(ray.org);
      const sse3f D = sse3f(ray.dir);
      const sse3f C = sse3f(tri.v0) - O;
      const sse3f R = cross(D,C);
      const ssef det = dot(sse3f(tri.Ng),D);
      const ssef absDet = abs(det);
      const ssef sgnDet = signmsk(det);

      /* perform edge tests */
      const ssef U = dot(R,sse3f(tri.e2)) ^ sgnDet;
      const ssef V = dot(R,sse3f(tri.e1)) ^ sgnDet;
      sseb valid = (det != ssef(zero)) & (U >= 0.0f) & (V >= 0.0f) & (U+V<=absDet);
      if (likely(none(valid))) return;
      
      /* perform depth test */
      const ssef T = dot(sse3f(tri.Ng),C) ^ sgnDet;
      valid &= (T > absDet*ssef(ray.tnear)) & (T < absDet*ssef(ray.tfar));
      if (likely(none(valid))) return;

      /* update hit information */
      const ssef rcpAbsDet = rcp(absDet);
      const ssef u = U * rcpAbsDet;
      const ssef v = V * rcpAbsDet;
      const ssef t = T * rcpAbsDet;
      const size_t i = select_min(valid,t);
      ray.u   = u[i];
      ray.v   = v[i];
      ray.tfar = t[i];
      ray.Ng.x = tri.Ng.x[i];
      ray.Ng.y = tri.Ng.y[i];
      ray.Ng.z = tri.Ng.z[i];
      ray.id0 = tri.id0[i];
      ray.id1 = tri.id1[i];
    }

    /*! Test if the ray is occluded by one of the triangles. */
    static __forceinline bool occluded(const Ray& ray, const Triangle4& tri, const Vec3fa* vertices = NULL)
    {
      STAT3(shadow.trav_tris,1,1,1);

      /* calculate determinant */
      const sse3f O = sse3f(ray.org);
      const sse3f D = sse3f(ray.dir);
      const sse3f C = sse3f(tri.v0) - O;
      const sse3f R = cross(D,C);
      const ssef det = dot(sse3f(tri.Ng),D);
      const ssef absDet = abs(det);
      const ssef sgnDet = signmsk(det);

      /* perform edge tests */
      const ssef U = dot(R,sse3f(tri.e2)) ^ sgnDet;
      const ssef V = dot(R,sse3f(tri.e1)) ^ sgnDet;
      const ssef W = absDet-U-V;
      sseb valid = (U >= 0.0f) & (V >= 0.0f) & (W >= 0.0f);
      if (unlikely(none(valid))) return false;
      
      /* perform depth test */
      const ssef T = dot(sse3f(tri.Ng),C) ^ sgnDet;
      valid &= (det != ssef(zero)) & (T >= absDet*ssef(ray.tnear)) & (absDet*ssef(ray.tfar) >= T);
      if (unlikely(none(valid))) return false;

      return true;
    }

    /*! Intersect a ray with the 4 triangles and updates the hit. */
    template<typename Ray>
      static __forceinline void intersect(const size_t k, Ray& ray, const Triangle4& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,1,1);

      /* calculate determinant */
      const sse3f O = sse3f(ray.org.x[k],ray.org.y[k],ray.org.z[k]);
      const sse3f D = sse3f(ray.dir.x[k],ray.dir.y[k],ray.dir.z[k]);
      const sse3f C = sse3f(tri.v0) - O;
      const sse3f R = cross(D,C);
      const ssef det = dot(sse3f(tri.Ng),D);
      const ssef absDet = abs(det);
      const ssef sgnDet = signmsk(det);

      /* perform edge tests */
      const ssef U = dot(R,sse3f(tri.e2)) ^ sgnDet;
      const ssef V = dot(R,sse3f(tri.e1)) ^ sgnDet;
      sseb valid = (det != ssef(zero)) & (U >= 0.0f) & (V >= 0.0f) & (U+V<=absDet);
      if (unlikely(none(valid))) return;
      
      /* perform depth test */
      const ssef T = dot(sse3f(tri.Ng),C) ^ sgnDet;
      valid &= (T > absDet*ssef(ray.tnear[k])) & (T < absDet*ssef(ray.tfar[k]));
      if (unlikely(none(valid))) return;

      /* update hit information */
      const ssef rcpAbsDet = rcp(absDet);
      const ssef u = U * rcpAbsDet;
      const ssef v = V * rcpAbsDet;
      const ssef t = T * rcpAbsDet;
      const size_t i = select_min(valid,t);
      ray.u[k]   = u[i];
      ray.v[k]   = v[i];
      ray.tfar[k] = t[i];
      ray.id0[k] = tri.id0[i];
      ray.id1[k] = tri.id1[i];
      ray.Ng.x[k] = tri.Ng.x[i];
      ray.Ng.y[k] = tri.Ng.y[i];
      ray.Ng.z[k] = tri.Ng.z[i];
    }

    /*! Test if the ray is occluded by one of the triangles. */
    template<typename Ray>
      static __forceinline bool occluded(const size_t k, const Ray& ray, const Triangle4& tri, const Vec3fa* vertices = NULL)
    {
      STAT3(shadow.trav_tris,1,1,1);

      /* calculate determinant */
      const sse3f O = sse3f(ray.org.x[k],ray.org.y[k],ray.org.z[k]);
      const sse3f D = sse3f(ray.dir.x[k],ray.dir.y[k],ray.dir.z[k]);
      const sse3f C = sse3f(tri.v0) - O;
      const sse3f R = cross(D,C);
      const ssef det = dot(sse3f(tri.Ng),D);
      const ssef absDet = abs(det);
      const ssef sgnDet = signmsk(det);

      /* perform edge tests */
      const ssef U = dot(R,sse3f(tri.e2)) ^ sgnDet;
      const ssef V = dot(R,sse3f(tri.e1)) ^ sgnDet;
      const ssef W = absDet-U-V;
      sseb valid = (U >= 0.0f) & (V >= 0.0f) & (W >= 0.0f);
      if (unlikely(none(valid))) return false;
      
      /* perform depth test */
      const ssef T = dot(sse3f(tri.Ng),C) ^ sgnDet;
      valid &= (det != ssef(zero)) & (T >= absDet*ssef(ray.tnear[k])) & (absDet*ssef(ray.tfar[k]) >= T);
      if (unlikely(none(valid))) return false;

      return true;
    }
  };
}

#endif


