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

#ifndef __EMBREE_ACCEL_TRIANGLE1V_INTERSECTOR_MOELLER_H__
#define __EMBREE_ACCEL_TRIANGLE1V_INTERSECTOR_MOELLER_H__

#include "triangle1.h"
#include "../common/ray.h"
#include "../common/hit.h"

namespace embree
{
  /*! Moeller Trumbore ray/triangle intersector for a single ray with
   *  individual pre-gathered triangles. */
  struct Triangle1vIntersectorMoellerTrumbore
  {
    typedef Triangle1v Triangle;

    /*! Intersect a ray with the 4 triangles and updates the hit. */
    static __forceinline void intersect(const Ray& ray, Hit& hit, const Triangle1v& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,1,1);

      /* calculate edges, geometry normal, and determinant */
      const Vec3f p0 = tri.v0;
      const Vec3f p1 = tri.v1;
      const Vec3f p2 = tri.v2;
      const Vec3f e1 = p0-p1;
      const Vec3f e2 = p2-p0;
      const Vec3f Ng = cross(e1,e2);

      /* calculate determinant */
      const Vec3f C = p0 - ray.org;
      const Vec3f R = cross(ray.dir,C);
      const float det = dot(Ng,ray.dir);
      const float absDet = abs(det);
      const float sgnDet = signmsk(det);
      if (unlikely(det == 0.0f)) return;
      
      /* test against edge p2 p0 */
      const float U = xorf(dot(R,e2),sgnDet);
      if (likely(U < 0.0f)) return;

      /* test against edge p0 p1 */
      const float V = xorf(dot(R,e1),sgnDet);
      if (likely(V < 0.0f)) return;

      /* test against edge p1 p2 */
      const float W = absDet-U-V;
      if (likely(W < 0.0f)) return;

      /* perform depth test */
      const float T = xorf(dot(Ng,C),sgnDet);
      if (unlikely(T < absDet*ray.near || absDet*hit.t < T)) return;

      /* update hit information */
      const float rcpAbsDet = rcp(absDet);
      hit.u   = U * rcpAbsDet;
      hit.v   = V * rcpAbsDet;
      hit.t   = T * rcpAbsDet;
      hit.id0 = tri.v0.a;
      hit.id1 = tri.v1.a;
    }

    /*! Test if the ray is occluded by one of the triangles. */
    static __forceinline bool occluded(const Ray& ray, const Triangle1v& tri, const Vec3fa* vertices = NULL)
    {
      STAT3(shadow.trav_tris,1,1,1);

      /* calculate edges, geometry normal, and determinant */
      const Vec3f p0 = tri.v0;
      const Vec3f p1 = tri.v1;
      const Vec3f p2 = tri.v2;
      const Vec3f e1 = p0-p1;
      const Vec3f e2 = p2-p0;
      const Vec3f Ng = cross(e1,e2);
      
      /* calculate determinant */
      const Vec3f C = p0 - ray.org;
      const Vec3f R = cross(ray.dir,C);
      const float det = dot(Ng,ray.dir);
      const float absDet = abs(det);
      const float sgnDet = signmsk(det);
      if (unlikely(det == 0.0f)) return false;

      /* test against edge p2 p0 */
      const float U = xorf(dot(R,e2),sgnDet);
      if (likely(U < 0.0f)) return false;

      /* test against edge p0 p1 */
      const float V = xorf(dot(R,e1),sgnDet);
      if (likely(V < 0.0f)) return false;

      /* test against edge p1 p2 */
      const float W = absDet-U-V;
      if (likely(W < 0.0f)) return false;

      /* perform depth test */
      const float T = xorf(dot(Ng,C),sgnDet);
      if (unlikely(T < absDet*ray.near || absDet*ray.far < T)) return false;

      return true;
    }
  };
}

#endif


