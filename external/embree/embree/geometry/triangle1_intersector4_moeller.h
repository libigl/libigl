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

#ifndef __EMBREE_ACCEL_TRIANGLE1_INTERSECTOR4_MOELLER_H__
#define __EMBREE_ACCEL_TRIANGLE1_INTERSECTOR4_MOELLER_H__

#include "triangle1.h"
#include "../common/ray4.h"

namespace embree
{
  /*! Intersector for individual precomputed triangles with 4
   *  rays. This intersector implements a modified version of the
   *  Moeller Trumbore intersector from the paper "Fast, Minimum
   *  Storage Ray-Triangle Intersection". In contrast to the paper we
   *  precalculate some factors and factor the calculations
   *  differently to allow precalculating the cross product e1 x
   *  e2. */
  struct Triangle1Intersector4MoellerTrumbore
  {
    typedef Triangle1 Triangle;
    
    /*! Name of intersector */
    static const char* name() { return "moeller"; }

    /*! Intersects a 4 rays with 4 triangles. */
    static __forceinline void intersect(const sseb& valid_i, Ray4& ray, const Triangle1& tri, const Vec3fa* vertices = NULL)
    {
      STAT3(normal.trav_tris,1,popcnt(valid_i),4);
      sseb valid = valid_i;
        
      /* calculate determinant */
      const sse3f C = sse3f(tri.v0) - ray.org;
      const sse3f R = cross(ray.dir,C);
      const ssef det = dot(sse3f(tri.Ng),ray.dir);
      const ssef absDet = abs(det);
      const ssef sgnDet = signmsk(det);
      valid &= det != ssef(zero);
      if (unlikely(none(valid))) return;
      
      /* test against edge p2 p0 */
      const ssef U = dot(R,sse3f(tri.e2)) ^ sgnDet;
      valid &= U >= 0.0f;
      if (likely(none(valid))) return;
      
      /* test against edge p0 p1 */
      const ssef V = dot(R,sse3f(tri.e1)) ^ sgnDet;
      valid &= V >= 0.0f;
      if (likely(none(valid))) return;
      
      /* test against edge p1 p2 */
      const ssef W = absDet-U-V;
      valid &= W >= 0.0f;
      if (likely(none(valid))) return;
      
      /* perform depth test */
      const ssef T = dot(sse3f(tri.Ng),C) ^ sgnDet;
      valid &= (T >= absDet*ray.tnear) & (absDet*ray.tfar >= T);
      if (unlikely(none(valid))) return;
      
      /* update hit information for all rays that hit the triangle */
      const ssef rcpAbsDet = rcp(absDet);
      ray.u   = select(valid,U * rcpAbsDet,ray.u );
      ray.v   = select(valid,V * rcpAbsDet,ray.v );
      ray.tfar = select(valid,T * rcpAbsDet,ray.tfar );
      ray.id0 = select(valid,tri.e1.a,ray.id0);
      ray.id1 = select(valid,tri.e2.a,ray.id1);
      ray.Ng.x = select(valid,tri.Ng.x,ray.Ng.x);
      ray.Ng.y = select(valid,tri.Ng.y,ray.Ng.y);
      ray.Ng.z = select(valid,tri.Ng.z,ray.Ng.z);
    }

    /*! Test for 4 rays if they are occluded by any of the 4 triangle. */
    static __forceinline sseb occluded(const sseb& valid_i, const Ray4& ray, const Triangle1& tri, const Vec3fa* vertices = NULL)
    {
      STAT3(shadow.trav_tris,1,popcnt(valid_i),4);

      /* load edges and geometry normal */
      sseb valid = valid_i;
      const Vector3f p0(tri.v0.x,tri.v0.y,tri.v0.z);
      const Vector3f e1(tri.e1.x,tri.e1.y,tri.e1.z);
      const Vector3f e2(tri.e2.x,tri.e2.y,tri.e2.z);
      const Vector3f Ng(tri.Ng.x,tri.Ng.y,tri.Ng.z);
        
      /* calculate determinant */
      const sse3f C = sse3f(p0) - ray.org;
      const sse3f R = cross(ray.dir,C);
      const ssef det = dot(sse3f(Ng),ray.dir);
      const ssef absDet = abs(det);
      const ssef sgnDet = signmsk(det);
      valid &= det != ssef(zero);
      if (unlikely(none(valid))) return valid;
      
      /* test against edge p2 p0 */
      const ssef U = dot(R,sse3f(e2)) ^ sgnDet;
      valid &= U >= 0.0f;
      if (likely(none(valid))) return valid;
      
      /* test against edge p0 p1 */
      const ssef V = dot(R,sse3f(e1)) ^ sgnDet;
      valid &= V >= 0.0f;
      if (likely(none(valid))) return valid;
      
      /* test against edge p1 p2 */
      const ssef W = absDet-U-V;
      valid &= W >= 0.0f;
      if (likely(none(valid))) return valid;
      
      /* perform depth test */
      const ssef T = dot(sse3f(Ng),C) ^ sgnDet;
      valid &= (T >= absDet*ray.tnear) & (absDet*ray.tfar >= T);
      return valid;
    }
  };
}

#endif


