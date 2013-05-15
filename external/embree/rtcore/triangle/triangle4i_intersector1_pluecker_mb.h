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

#ifndef __EMBREE_ACCEL_TRIANGLE4I_INTERSECTOR_PLUECKER_MB_H__
#define __EMBREE_ACCEL_TRIANGLE4I_INTERSECTOR_PLUECKER_MB_H__

#include "triangle4i.h"
#include "../common/ray.h"
#include "../common/hit.h"

namespace embree
{
  /*! Pluecker intersector for 4 triangles from an indexed face set
   *  with a single ray. The test first shifts the ray origin into the
   *  origin of the coordinate system and then uses Pluecker
   *  coordinates for the intersection. Due to the shift, the Pluecker
   *  coordinate calculation simplifies. */
  struct Triangle4iIntersectorPlueckerMB
  {
    typedef Triangle4i Triangle;

    /*! Intersect a ray with the 4 triangles and updates the hit. */
    static __forceinline void intersect(const Ray& ray, Hit& hit, const Triangle4i& tri, const Vec3fa* vertices = NULL)
    {
      STAT3(normal.trav_tris,1,1,1);

      /* gather vertices */
      const ssef* vert = (const ssef*) vertices;
      sse3f p0; transpose(vert[tri.v0[0]],vert[tri.v0[1]],vert[tri.v0[2]],vert[tri.v0[3]], p0.x,p0.y,p0.z);
      sse3f p1; transpose(vert[tri.v1[0]],vert[tri.v1[1]],vert[tri.v1[2]],vert[tri.v1[3]], p1.x,p1.y,p1.z);
      sse3f p2; transpose(vert[tri.v2[0]],vert[tri.v2[1]],vert[tri.v2[2]],vert[tri.v2[3]], p2.x,p2.y,p2.z);
      
      const sseb motion = (tri.id0 & ssei(0x80000000)) != ssei(zero);
      if (unlikely(any(tri.valid() & motion))) {
        const ssei dv0 = select(motion,tri.v0+1,tri.v0);
        const ssei dv1 = select(motion,tri.v1+1,tri.v1);
        const ssei dv2 = select(motion,tri.v2+1,tri.v2);
        sse3f dp0; transpose(vert[dv0[0]],vert[dv0[1]],vert[dv0[2]],vert[dv0[3]], dp0.x,dp0.y,dp0.z);
        sse3f dp1; transpose(vert[dv1[0]],vert[dv1[1]],vert[dv1[2]],vert[dv1[3]], dp1.x,dp1.y,dp1.z);
        sse3f dp2; transpose(vert[dv2[0]],vert[dv2[1]],vert[dv2[2]],vert[dv2[3]], dp2.x,dp2.y,dp2.z);
        p0 = select(motion,p0+ssef(ray.time)*dp0,p0);
        p1 = select(motion,p1+ssef(ray.time)*dp1,p1);
        p2 = select(motion,p2+ssef(ray.time)*dp2,p2);
      }
      
      /* calculate edges, geometry normal, and determinant */
      const sse3f O = ray.org;
      const sse3f D = ray.dir;
      const sse3f v0 = p0-O;
      const sse3f v1 = p1-O;
      const sse3f v2 = p2-O;
      const sse3f e0 = v2-v0;
      const sse3f e1 = v0-v1;
      const sse3f e2 = v1-v2;
      const sse3f Ng = ssef(2.0f)*cross(e1,e0);
      const ssef det = dot(Ng,D);
      const ssef absDet = abs(det);
      const ssef sgnDet = signmsk(det);

      /* perform edge tests */
      const ssef U = dot(cross(v2+v0,e0),D) ^ sgnDet;
      const ssef V = dot(cross(v0+v1,e1),D) ^ sgnDet;
      const ssef W = dot(cross(v1+v2,e2),D) ^ sgnDet;
      sseb valid = tri.valid() & (U >= 0.0f) & (V >= 0.0f) & (W >= 0.0f);
      if (unlikely(none(valid))) return;

      /* perform depth test */
      const ssef T = dot(v0,Ng) ^ sgnDet;
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
    static __forceinline bool occluded(const Ray& ray, const Triangle4i& tri, const Vec3fa* vertices = NULL)
    {
      STAT3(shadow.trav_tris,1,1,1);

      /* gather vertices */
      const ssef* vert = (const ssef*) vertices;
      sse3f p0; transpose(vert[tri.v0[0]],vert[tri.v0[1]],vert[tri.v0[2]],vert[tri.v0[3]], p0.x,p0.y,p0.z);
      sse3f p1; transpose(vert[tri.v1[0]],vert[tri.v1[1]],vert[tri.v1[2]],vert[tri.v1[3]], p1.x,p1.y,p1.z);
      sse3f p2; transpose(vert[tri.v2[0]],vert[tri.v2[1]],vert[tri.v2[2]],vert[tri.v2[3]], p2.x,p2.y,p2.z);
      
      const sseb motion = (tri.id0 & ssei(0x80000000)) != ssei(zero);
      if (unlikely(any(tri.valid() & motion))) {
        const ssei dv0 = select(motion,tri.v0+1,tri.v0);
        const ssei dv1 = select(motion,tri.v1+1,tri.v1);
        const ssei dv2 = select(motion,tri.v2+1,tri.v2);
        sse3f dp0; transpose(vert[dv0[0]],vert[dv0[1]],vert[dv0[2]],vert[dv0[3]], dp0.x,dp0.y,dp0.z);
        sse3f dp1; transpose(vert[dv1[0]],vert[dv1[1]],vert[dv1[2]],vert[dv1[3]], dp1.x,dp1.y,dp1.z);
        sse3f dp2; transpose(vert[dv2[0]],vert[dv2[1]],vert[dv2[2]],vert[dv2[3]], dp2.x,dp2.y,dp2.z);
        p0 = select(motion,p0+ssef(ray.time)*dp0,p0);
        p1 = select(motion,p1+ssef(ray.time)*dp1,p1);
        p2 = select(motion,p2+ssef(ray.time)*dp2,p2);
      }
      
      /* calculate edges, geometry normal, and determinant */
      const sse3f O = ray.org;
      const sse3f D = ray.dir;
      const sse3f v0 = p0-O;
      const sse3f v1 = p1-O;
      const sse3f v2 = p2-O;
      const sse3f e0 = v2-v0;
      const sse3f e1 = v0-v1;
      const sse3f e2 = v1-v2;
      const sse3f Ng = ssef(2.0f)*cross(e1,e0);
      const ssef det = dot(Ng,D);
      const ssef absDet = abs(det);
      const ssef sgnDet = signmsk(det);

      /* perform edge tests */
      const ssef U = dot(cross(v2+v0,e0),D) ^ sgnDet;
      const ssef V = dot(cross(v0+v1,e1),D) ^ sgnDet;
      const ssef W = dot(cross(v1+v2,e2),D) ^ sgnDet;
      sseb valid = tri.valid() & (U >= 0.0f) & (V >= 0.0f) & (W >= 0.0f);
      if (unlikely(none(valid))) return false;
      
      /* perform depth test */
      const ssef T = dot(v0,Ng) ^ sgnDet;
      valid &= (det != ssef(zero)) & (T >= absDet*ssef(ray.near)) & (absDet*ssef(ray.far) >= T);
      if (unlikely(none(valid))) return false;

      return true;
    }
  };
}

#endif
