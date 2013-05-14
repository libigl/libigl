// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#ifndef __EMBREE_ACCEL_TRIANGLE8_INTERSECTOR1_MOELLER_H__
#define __EMBREE_ACCEL_TRIANGLE8_INTERSECTOR1_MOELLER_H__

#include "triangle8.h"
#include "../common/ray.h"
#include "../common/hit.h"

namespace embree
{
  /*! Intersector for a single ray with 8 triangles. This intersector
   *  implements a modified version of the Moeller Trumbore
   *  intersector from the paper "Fast, Minimum Storage Ray-Triangle
   *  Intersection". In contrast to the paper we precalculate some
   *  factors and factor the calculations differently to allow
   *  precalculating the cross product e1 x e2. */
  struct Triangle8IntersectorMoellerTrumbore
  {
    typedef Triangle8 Triangle;

    /*! Intersect a ray with the 4 triangles and updates the hit. */
    static __forceinline void intersect(const Ray& ray, Hit& hit, const Triangle8& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,1,1);

      /* calculate determinant */
      const avx3f O = avx3f(ray.org);
      const avx3f D = avx3f(ray.dir);
      const avx3f C = tri.v0 - O;
      const avx3f R = cross(D,C);
      const avxf det = dot(tri.Ng,D);
      const avxf absDet = abs(det);
      const avxf sgnDet = signmsk(det);

      /* perform edge tests */
      const avxf U = dot(R,tri.e2) ^ sgnDet;
      const avxf V = dot(R,tri.e1) ^ sgnDet;
      avxb valid = (det != avxf(zero)) & (U >= 0.0f) & (V >= 0.0f) & (U+V<=absDet);
      if (unlikely(none(valid))) return;
      
      /* perform depth test */
      const avxf T = dot(tri.Ng,C) ^ sgnDet;
      valid &= (T > absDet*avxf(ray.near)) & (T < absDet*avxf(hit.t));
      if (unlikely(none(valid))) return;

      /* update hit information */
      const avxf rcpAbsDet = rcp(absDet);
      const avxf u = U * rcpAbsDet;
      const avxf v = V * rcpAbsDet;
      const avxf t = T * rcpAbsDet;
      const size_t i = select_min(valid,t);
      hit.u   = u[i];
      hit.v   = v[i];
      hit.t   = t[i];
      hit.id0 = tri.id0[i];
      hit.id1 = tri.id1[i];
    }

    /*! Test if the ray is occluded by one of the triangles. */
    static __forceinline bool occluded(const Ray& ray, const Triangle8& tri, const Vec3fa* vertices = NULL)
    {
      STAT3(shadow.trav_tris,1,1,1);

      /* calculate determinant */
      const avx3f O = avx3f(ray.org);
      const avx3f D = avx3f(ray.dir);
      const avx3f C = tri.v0 - O;
      const avx3f R = cross(D,C);
      const avxf det = dot(tri.Ng,D);
      const avxf absDet = abs(det);
      const avxf sgnDet = signmsk(det);

      /* perform edge tests */
      const avxf U = dot(R,tri.e2) ^ sgnDet;
      const avxf V = dot(R,tri.e1) ^ sgnDet;
      const avxf W = absDet-U-V;
      avxb valid = (U >= 0.0f) & (V >= 0.0f) & (W >= 0.0f);
      if (unlikely(none(valid))) return false;
      
      /* perform depth test */
      const avxf T = dot(tri.Ng,C) ^ sgnDet;
      valid &= (det != avxf(zero)) & (T >= absDet*avxf(ray.near)) & (absDet*avxf(ray.far) >= T);
      if (unlikely(none(valid))) return false;

      return true;
    }
  };
}

#endif


