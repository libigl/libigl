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

#ifndef __EMBREE_ACCEL_TRIANGLE4XFM_INTERSECTOR_MOELLER_H__
#define __EMBREE_ACCEL_TRIANGLE4XFM_INTERSECTOR_MOELLER_H__

#include "triangle4xfm.h"
#include "../common/ray.h"
#include "../common/hit.h"

namespace embree
{
  /*! Intersector for a single ray with 4 triangles. Performs the
   *  ray/triangle intersection in a coordinate space where the
   *  triangle is a unit triangle. */
  struct Triangle4XfmIntersectorMoellerTrumbore
  {
    typedef Triangle4Xfm Triangle;

    /*! Intersect a ray with the 4 triangles and updates the hit. */
    static __forceinline void intersect(const Ray& ray, Hit& hit, const Triangle4Xfm& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,1,1);
      const sse3f org = ray.org;
      const sse3f dir = ray.dir;
      const Triangle4Xfm::AffineSpace4x3f& xfm = tri.xfm;
      
      /* transform ray into unit triangle space */
      const ssef uorgz = org.x*xfm.l.vx.z + org.y*xfm.l.vy.z + org.z*xfm.l.vz.z + xfm.p.z;
      const ssef udirz = dir.x*xfm.l.vx.z + dir.y*xfm.l.vy.z + dir.z*xfm.l.vz.z;
      const ssef uorgy = org.x*xfm.l.vx.y + org.y*xfm.l.vy.y + org.z*xfm.l.vz.y + xfm.p.y;
      const ssef udiry = dir.x*xfm.l.vx.y + dir.y*xfm.l.vy.y + dir.z*xfm.l.vz.y;
      const ssef uorgx = org.x*xfm.l.vx.x + org.y*xfm.l.vy.x + org.z*xfm.l.vz.x + xfm.p.x;
      const ssef udirx = dir.x*xfm.l.vx.x + dir.y*xfm.l.vy.x + dir.z*xfm.l.vz.x;
      const ssef absudirz = abs(udirz), sgnudirz = signmsk(udirz);

      /* perform edge tests */
      const ssef V = (uorgy*udirz-uorgz*udiry)^sgnudirz;
      const ssef U = (uorgx*udirz-uorgz*udirx)^sgnudirz;
      sseb valid = (V >= 0.0f) & (U >= 0.0f) & (U+V <= absudirz);
      if (unlikely(none(valid))) return;
      
      /* perform depth test */
      const ssef T = -uorgz^sgnudirz;
      valid &= (udirz != 0.0f) & (T > absudirz*ssef(ray.near)) & (T < absudirz*ssef(hit.t));
      if (unlikely(none(valid))) return;

      /* update hit information */
      const ssef rcpabsudirz = rcp(absudirz);
      const ssef u = U * rcpabsudirz;
      const ssef v = V * rcpabsudirz;
      const ssef t = T * rcpabsudirz;
      const size_t i = select_min(valid,t);
      hit.u   = u[i];
      hit.v   = v[i];
      hit.t   = t[i];
      hit.id0 = tri.id0[i];
      hit.id1 = tri.id1[i];
    }

    /*! Test if the ray is occluded by one of the triangles. */
    static __forceinline bool occluded(const Ray& ray, const Triangle4Xfm& tri, const Vec3fa* vertices = NULL)
    {
      STAT3(shadow.trav_tris,1,1,1);
      const sse3f org = ray.org;
      const sse3f dir = ray.dir;
      const Triangle4Xfm::AffineSpace4x3f& xfm = tri.xfm;
      
      /* transform ray into unit triangle space */
      const ssef uorgz = org.x*xfm.l.vx.z + org.y*xfm.l.vy.z + org.z*xfm.l.vz.z + xfm.p.z;
      const ssef udirz = dir.x*xfm.l.vx.z + dir.y*xfm.l.vy.z + dir.z*xfm.l.vz.z;
      const ssef uorgy = org.x*xfm.l.vx.y + org.y*xfm.l.vy.y + org.z*xfm.l.vz.y + xfm.p.y;
      const ssef udiry = dir.x*xfm.l.vx.y + dir.y*xfm.l.vy.y + dir.z*xfm.l.vz.y;
      const ssef uorgx = org.x*xfm.l.vx.x + org.y*xfm.l.vy.x + org.z*xfm.l.vz.x + xfm.p.x;
      const ssef udirx = dir.x*xfm.l.vx.x + dir.y*xfm.l.vy.x + dir.z*xfm.l.vz.x;
      const ssef absudirz = abs(udirz), sgnudirz = signmsk(udirz);

      /* perform edge tests */
      const ssef V = (uorgy*udirz-uorgz*udiry)^sgnudirz;
      const ssef U = (uorgx*udirz-uorgz*udirx)^sgnudirz;
      sseb valid = (V >= 0.0f) & (U >= 0.0f) & (U+V <= absudirz);
      if (unlikely(none(valid))) return false;
      
      /* perform depth test */
      const ssef T = -uorgz^sgnudirz;
      valid &= (udirz != 0.0f) & (T > absudirz*ssef(ray.near)) & (T < absudirz*ssef(ray.far));
      return any(valid);
    }
  };
}

#endif


