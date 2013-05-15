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

#ifndef __EMBREE_ACCEL_TRIANGLE1XFM_INTERSECTOR1_MOELLER_H__
#define __EMBREE_ACCEL_TRIANGLE1XFM_INTERSECTOR1_MOELLER_H__

#include "triangle1xfm.h"
#include "../common/ray.h"
#include "../common/hit.h"

namespace embree
{
  /*! Intersector for a single ray with a single triangle. This
   *  intersector performs the ray/triangle intersection in a
   *  coordinate space where the triangle is a unit triangle as
   *  described in Woop's master thesis "A Ray Tracing Hardware
   *  Architecture for Dynamic Scenes". */
  struct Triangle1XfmIntersectorMoellerTrumbore
  {
    typedef Triangle1Xfm Triangle;

    /*! Intersect a ray with the 4 triangles and updates the hit. */
    static __forceinline void intersect(const Ray& ray, Hit& hit, const Triangle1Xfm& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,1,1);
      
      /* transform ray into unit triangle space */
      const Vec3f uorg = xfmPoint (tri.xfm,ray.org);
      const Vec3f udir = xfmVector(tri.xfm,ray.dir);
      const float absudirz = abs(udir.z), sgnudirz = signmsk(udir.z);

      /* perform edge tests */
      const float V = xorf(uorg.y*udir.z-uorg.z*udir.y,sgnudirz);
      const float U = xorf(uorg.x*udir.z-uorg.z*udir.x,sgnudirz);
      sseb valid = (V >= 0.0f) & (U >= 0.0f) & (U+V <= absudirz);
      if (unlikely(none(valid))) return;
      
      /* perform depth test */
      const float T = -xorf(uorg.z,sgnudirz);
      valid &= (udir.z != 0.0f) & (T > absudirz*ray.near) & (T < absudirz*hit.t);
      if (unlikely(none(valid))) return;

      /* update hit information */
      const float rcpabsudirz = rcp(absudirz);
      hit.u   = U * rcpabsudirz;
      hit.v   = V * rcpabsudirz;
      hit.t   = T * rcpabsudirz;
      hit.id0 = tri.id0;
      hit.id1 = tri.id1;
    }

    /*! Test if the ray is occluded by one of the triangles. */
    static __forceinline bool occluded(const Ray& ray, const Triangle1Xfm& tri, const Vec3fa* vertices = NULL)
    {
      STAT3(shadow.trav_tris,1,1,1);
      
      /* transform ray into unit triangle space */
      const Vec3f uorg = xfmPoint (tri.xfm,ray.org);
      const Vec3f udir = xfmVector(tri.xfm,ray.dir);
      const float absudirz = abs(udir.z), sgnudirz = signmsk(udir.z);
      
      /* perform edge tests */
      const float V = xorf(uorg.y*udir.z-uorg.z*udir.y,sgnudirz);
      const float U = xorf(uorg.x*udir.z-uorg.z*udir.x,sgnudirz);
      sseb valid = (V >= 0.0f) & (U >= 0.0f) & (U+V <= absudirz);
      if (unlikely(none(valid))) return false;
      
      /* perform depth test */
      const float T = -xorf(uorg.z,sgnudirz);
      valid &= (udir.z != 0.0f) & (T > absudirz*ray.near) & (T < absudirz*ray.far);
      return any(valid);
    }
  };
}

#endif


