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

#ifndef __EMBREE_ACCEL_TRIANGLE4V_INTERSECTOR_MOELLER_H__
#define __EMBREE_ACCEL_TRIANGLE4V_INTERSECTOR_MOELLER_H__

#include "triangle4v.h"
#include "../common/ray.h"
#include "../common/hit.h"

namespace embree
{
  /*! Intersector for 4 triangles pre-gathered triangles with a single
   *  ray. Implements a modified version of the Moeller Trumbore
   *  intersector. */
  struct Triangle4vIntersectorMoellerTrumbore
  {
    typedef Triangle4v Triangle;

    /*! Intersect a ray with the 4 triangles and updates the hit. */
    static __forceinline void intersect(const Ray& ray, Hit& hit, const Triangle4v& tri, const Vec3fa* vertices = NULL)
    {
      STAT3(normal.trav_tris,1,1,1);

      /* calculate edges and geometry normal */
      const sse3f p0 = tri.v0;
      const sse3f p1 = tri.v1;
      const sse3f p2 = tri.v2;
      const sse3f e1 = p0-p1;
      const sse3f e2 = p2-p0;
      const sse3f Ng = cross(e1,e2);
      const sse3f O = sse3f(ray.org);
      const sse3f D = sse3f(ray.dir);
      
      /* calculate determinant */
      const sse3f C = p0 - O;
      const sse3f R = cross(D,C);
      const ssef det = dot(Ng,D);
      const ssef absDet = abs(det);
      const ssef sgnDet = signmsk(det);

      /* perform edge tests */
      const ssef U = dot(R,e2) ^ sgnDet;
      const ssef V = dot(R,e1) ^ sgnDet;
      const ssef W = absDet-U-V;
      sseb valid = tri.valid() & (U >= 0.0f) & (V >= 0.0f) & (W >= 0.0f);
      if (unlikely(none(valid))) return;
      
      /* perform depth test */
      const ssef T = dot(Ng,C) ^ sgnDet;
      valid &= (det != ssef(zero)) & (T >= absDet*ssef(ray.near)) & (absDet*ssef(hit.t) >= T);
      if (unlikely(none(valid))) return;

      /* update hit information */
      const ssef rcpAbsDet = rcp(absDet);
      const ssef u = U * rcpAbsDet;
      const ssef v = V * rcpAbsDet;
      const ssef t = T * rcpAbsDet;
      const size_t i = select_min(valid,t);
      hit.t = t[i];
      hit.u = u[i];
      hit.v = v[i];
      hit.id0 = tri.id0[i] & 0x7FFFFFFF;
      hit.id1 = tri.id1[i];
    }

    /*! Test if the ray is occluded by one of the triangles. */
    static __forceinline bool occluded(const Ray& ray, const Triangle4v& tri, const Vec3fa* vertices = NULL)
    {
      STAT3(shadow.trav_tris,1,1,1);

      /* calculate edges and geometry normal */
      const sse3f p0 = tri.v0;
      const sse3f p1 = tri.v1;
      const sse3f p2 = tri.v2;
      const sse3f e1 = p0-p1;
      const sse3f e2 = p2-p0;
      const sse3f Ng = cross(e1,e2);
      const sse3f O = sse3f(ray.org);
      const sse3f D = sse3f(ray.dir);
      
      /* calculate determinant */
      const sse3f C = p0 - O;
      const sse3f R = cross(D,C);
      const ssef det = dot(Ng,D);
      const ssef absDet = abs(det);
      const ssef sgnDet = signmsk(det);

      /* perform edge tests */
      const ssef U = dot(R,e2) ^ sgnDet;
      const ssef V = dot(R,e1) ^ sgnDet;
      const ssef W = absDet-U-V;
      sseb valid = tri.valid() & (U >= 0.0f) & (V >= 0.0f) & (W >= 0.0f);
      if (unlikely(none(valid))) return false;
      
      /* perform depth test */
      const ssef T = dot(Ng,C) ^ sgnDet;
      valid &= (det != ssef(zero)) & (T >= absDet*ssef(ray.near)) & (absDet*ssef(ray.far) >= T);
      if (unlikely(none(valid))) return false;

      return true;
    }
  };
}

#endif
