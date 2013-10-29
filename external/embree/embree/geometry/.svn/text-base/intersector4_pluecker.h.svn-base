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

#ifndef __EMBREE_INTERSECTOR4_PLUECKER_H__
#define __EMBREE_INTERSECTOR4_PLUECKER_H__

#include "../common/ray4.h"

namespace embree
{
  /*! Modified Pluecker ray/triangle intersector. The test first shifts the ray
   *  origin into the origin of the coordinate system and then uses
   *  Pluecker coordinates for the intersection. Due to the shift, the
   *  Pluecker coordinate calculation simplifies. The edge equations
   *  are watertight along the edge for neighboring triangles. */
  struct Intersector4Pluecker
  {
    /*! Name of intersector */
    static const char* name() { return "pluecker"; }

    /*! Intersect a ray with the 4 triangles and updates the hit. */
    static __forceinline void intersect(const sseb& valid_i, Ray4& ray, const Vector3f& p0, const Vector3f& p1, const Vector3f& p2, const int& id0, const int& id1)
    {
      /* calculate vertices relative to ray origin */
      sseb valid = valid_i;
      const sse3f O = ray.org;
      const sse3f D = ray.dir;
      const sse3f v0 = sse3f(p0)-O;
      const sse3f v1 = sse3f(p1)-O;
      const sse3f v2 = sse3f(p2)-O;

      /* calculate triangle edges */
      const sse3f e0 = v2-v0;
      const sse3f e1 = v0-v1;
      const sse3f e2 = v1-v2;

      /* calculate geometry normal and determinant */
      const sse3f Ng = cross(e1,e0);
      const sse3f Ng2 = Ng+Ng;
      const ssef det = dot(sse3f(Ng2),D);
      const ssef absDet = abs(det);
      const ssef sgnDet = signmsk(det);
      valid &= det != ssef(zero);
      if (unlikely(none(valid))) return;

      /* perform edge tests */
      const ssef U = dot(sse3f(cross(v2+v0,e0)),D) ^ sgnDet;
      valid &= U >= 0.0f;
      if (likely(none(valid))) return;
      const ssef V = dot(sse3f(cross(v0+v1,e1)),D) ^ sgnDet;
      valid &= V >= 0.0f;
      if (likely(none(valid))) return;
      const ssef W = dot(sse3f(cross(v1+v2,e2)),D) ^ sgnDet;
      valid &= W >= 0.0f;
      if (likely(none(valid))) return;
      
      /* perform depth test */
      const ssef T = dot(v0,sse3f(Ng2)) ^ sgnDet;
      valid &= (T >= absDet*ray.tnear) & (absDet*ray.tfar >= T);
      if (unlikely(none(valid))) return;

      /* update hit information for all rays that hit the triangle */
      const ssef rcpAbsDet = rcp(absDet);
      ray.u   = select(valid,U * rcpAbsDet,ray.u );
      ray.v   = select(valid,V * rcpAbsDet,ray.v );
      ray.tfar = select(valid,T * rcpAbsDet,ray.tfar );
      ray.id0 = select(valid,id0,ray.id0);
      ray.id1 = select(valid,id1,ray.id1);
      ray.Ng.x = select(valid,Ng2.x,ray.Ng.x);
      ray.Ng.y = select(valid,Ng2.y,ray.Ng.y);
      ray.Ng.z = select(valid,Ng2.z,ray.Ng.z);
    }

    /*! Intersect 4 rays with 4 triangles and updates the hit. */
    static __forceinline void intersect(const sseb& valid_i, Ray4& ray, 
                                        const sse3f& p0, const sse3f& p1, const sse3f& p2, const ssei& id0, const ssei& id1)
    {
      /* calculate vertices relative to ray origin */
      sseb valid = valid_i;
      const sse3f O = ray.org;
      const sse3f D = ray.dir;
      const sse3f v0 = p0-O;
      const sse3f v1 = p1-O;
      const sse3f v2 = p2-O;

      /* calculate triangle edges */
      const sse3f e0 = v2-v0;
      const sse3f e1 = v0-v1;
      const sse3f e2 = v1-v2;

      /* calculate geometry normal and determinant */
      const sse3f Ng = cross(e1,e0);
      const sse3f Ng2 = Ng+Ng;
      const ssef det = dot(sse3f(Ng2),D);
      const ssef absDet = abs(det);
      const ssef sgnDet = signmsk(det);
      valid &= det != ssef(zero);
      if (unlikely(none(valid))) return;

      /* perform edge tests */
      const ssef U = dot(sse3f(cross(v2+v0,e0)),D) ^ sgnDet;
      valid &= U >= 0.0f;
      if (likely(none(valid))) return;
      const ssef V = dot(sse3f(cross(v0+v1,e1)),D) ^ sgnDet;
      valid &= V >= 0.0f;
      if (likely(none(valid))) return;
      const ssef W = dot(sse3f(cross(v1+v2,e2)),D) ^ sgnDet;
      valid &= W >= 0.0f;
      if (likely(none(valid))) return;
      
      /* perform depth test */
      const ssef T = dot(v0,sse3f(Ng2)) ^ sgnDet;
      valid &= (T >= absDet*ray.tnear) & (absDet*ray.tfar >= T);
      if (unlikely(none(valid))) return;

      /* update hit information for all rays that hit the triangle */
      const ssef rcpAbsDet = rcp(absDet);
      ray.u   = select(valid,U * rcpAbsDet,ray.u );
      ray.v   = select(valid,V * rcpAbsDet,ray.v );
      ray.tfar = select(valid,T * rcpAbsDet,ray.tfar );
      ray.id0 = select(valid,id0,ray.id0);
      ray.id1 = select(valid,id1,ray.id1);
      ray.Ng.x = select(valid,Ng2.x,ray.Ng.x);
      ray.Ng.y = select(valid,Ng2.y,ray.Ng.y);
      ray.Ng.z = select(valid,Ng2.z,ray.Ng.z);
    }

    /*! Test for all rays if they are occluded by the triangle. */
    static __forceinline sseb occluded(const sseb& valid_i, const Ray4& ray, const Vector3f& p0, const Vector3f& p1, const Vector3f& p2)
    {
      /* calculate vertices relative to ray origin */
      sseb valid = valid_i;
      const sse3f O = ray.org;
      const sse3f D = ray.dir;
      const sse3f v0 = sse3f(p0)-O;
      const sse3f v1 = sse3f(p1)-O;
      const sse3f v2 = sse3f(p2)-O;

      /* calculate triangle edges */
      const sse3f e0 = v2-v0;
      const sse3f e1 = v0-v1;
      const sse3f e2 = v1-v2;

      /* calculate geometry normal and determinant */
      const sse3f Ng = cross(e1,e0);
      const sse3f Ng2 = Ng+Ng;
      const ssef det = dot(sse3f(Ng2),D);
      const ssef absDet = abs(det);
      const ssef sgnDet = signmsk(det);
      valid &= det != ssef(zero);
      if (unlikely(none(valid))) return valid;

      /* perform edge tests */
      const ssef U = dot(sse3f(cross(v2+v0,e0)),D) ^ sgnDet;
      valid &= U >= 0.0f;
      if (likely(none(valid))) return valid;
      const ssef V = dot(sse3f(cross(v0+v1,e1)),D) ^ sgnDet;
      valid &= V >= 0.0f;
      if (likely(none(valid))) return valid;
      const ssef W = dot(sse3f(cross(v1+v2,e2)),D) ^ sgnDet;
      valid &= W >= 0.0f;
      if (likely(none(valid))) return valid;
      
      /* perform depth test */
      const ssef T = dot(v0,sse3f(Ng2)) ^ sgnDet;
      valid &= (T >= absDet*ray.tnear) & (absDet*ray.tfar >= T);
      return valid;
    }

    /*! Test for all rays if they are occluded by the triangle. */
    static __forceinline sseb occluded(const sseb& valid_i, const Ray4& ray, const sse3f& p0, const sse3f& p1, const sse3f& p2)
    {
      /* calculate vertices relative to ray origin */
      sseb valid = valid_i;
      const sse3f O = ray.org;
      const sse3f D = ray.dir;
      const sse3f v0 = p0-O;
      const sse3f v1 = p1-O;
      const sse3f v2 = p2-O;

      /* calculate triangle edges */
      const sse3f e0 = v2-v0;
      const sse3f e1 = v0-v1;
      const sse3f e2 = v1-v2;

      /* calculate geometry normal and determinant */
      const sse3f Ng = cross(e1,e0);
      const sse3f Ng2 = Ng+Ng;
      const ssef det = dot(sse3f(Ng2),D);
      const ssef absDet = abs(det);
      const ssef sgnDet = signmsk(det);
      valid &= det != ssef(zero);
      if (unlikely(none(valid))) return valid;

      /* perform edge tests */
      const ssef U = dot(sse3f(cross(v2+v0,e0)),D) ^ sgnDet;
      valid &= U >= 0.0f;
      if (likely(none(valid))) return valid;
      const ssef V = dot(sse3f(cross(v0+v1,e1)),D) ^ sgnDet;
      valid &= V >= 0.0f;
      if (likely(none(valid))) return valid;
      const ssef W = dot(sse3f(cross(v1+v2,e2)),D) ^ sgnDet;
      valid &= W >= 0.0f;
      if (likely(none(valid))) return valid;
      
      /* perform depth test */
      const ssef T = dot(v0,sse3f(Ng2)) ^ sgnDet;
      valid &= (T >= absDet*ray.tnear) & (absDet*ray.tfar >= T);
      return valid;
    }
  };
}

#endif


