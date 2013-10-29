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

#ifndef __EMBREE_INTERSECTOR8_PLUECKER_H__
#define __EMBREE_INTERSECTOR8_PLUECKER_H__

#include "../common/ray8.h"

namespace embree
{
  /*! Modified Pluecker ray/triangle intersector. The test first shifts the ray
   *  origin into the origin of the coordinate system and then uses
   *  Pluecker coordinates for the intersection. Due to the shift, the
   *  Pluecker coordinate calculation simplifies. The edge equations
   *  are watertight along the edge for neighboring triangles. */
  struct Intersector8Pluecker
  {
    /*! Name of intersector */
    static const char* name() { return "pluecker"; }

    /*! Intersect a ray with the 8 triangles and updates the hit. */
    static __forceinline void intersect(const avxb& valid_i, Ray8& ray,
                                        const Vector3f& p0, const Vector3f& p1, const Vector3f& p2, const int& id0, const int& id1)
    {
      /* calculate vertices relative to ray origin */
      avxb valid = valid_i;
      const avx3f O = ray.org;
      const avx3f D = ray.dir;
      const avx3f v0 = avx3f(p0)-O;
      const avx3f v1 = avx3f(p1)-O;
      const avx3f v2 = avx3f(p2)-O;

      /* calculate triangle edges */
      const avx3f e0 = v2-v0;
      const avx3f e1 = v0-v1;
      const avx3f e2 = v1-v2;

      /* calculate geometry normal and determinant */
      const avx3f Ng = cross(e1,e0);
      const avx3f Ng2 = Ng+Ng;
      const avxf det = dot(avx3f(Ng2),D);
      const avxf absDet = abs(det);
      const avxf sgnDet = signmsk(det);
      valid &= det != avxf(zero);
      if (unlikely(none(valid))) return;

      /* perform edge tests */
      const avxf U = dot(avx3f(cross(v2+v0,e0)),D) ^ sgnDet;
      valid &= U >= 0.0f;
      if (likely(none(valid))) return;
      const avxf V = dot(avx3f(cross(v0+v1,e1)),D) ^ sgnDet;
      valid &= V >= 0.0f;
      if (likely(none(valid))) return;
      const avxf W = dot(avx3f(cross(v1+v2,e2)),D) ^ sgnDet;
      valid &= W >= 0.0f;
      if (likely(none(valid))) return;
      
      /* perform depth test */
      const avxf T = dot(v0,avx3f(Ng2)) ^ sgnDet;
      valid &= (T >= absDet*ray.tnear) & (absDet*ray.tfar >= T);
      if (unlikely(none(valid))) return;

      /* update hit information for all rays that hit the triangle */
      const avxf rcpAbsDet = rcp(absDet);
      ray.u   = select(valid,U * rcpAbsDet,ray.u );
      ray.v   = select(valid,V * rcpAbsDet,ray.v );
      ray.tfar = select(valid,T * rcpAbsDet,ray.tfar );
      ray.id0 = select(valid,id0,ray.id0);
      ray.id1 = select(valid,id1,ray.id1);
      ray.Ng.x = select(valid,Ng2.x,ray.Ng.x);
      ray.Ng.y = select(valid,Ng2.y,ray.Ng.y);
      ray.Ng.z = select(valid,Ng2.z,ray.Ng.z);
    }

    /*! Intersect 8 rays with 8 triangles and updates the hit. */
    static __forceinline void intersect(const avxb& valid_i, Ray8& ray,
                                        const avx3f& p0, const avx3f& p1, const avx3f& p2, const avxi& id0, const avxi& id1)
    {
      /* calculate vertices relative to ray origin */
      avxb valid = valid_i;
      const avx3f O = ray.org;
      const avx3f D = ray.dir;
      const avx3f v0 = p0-O;
      const avx3f v1 = p1-O;
      const avx3f v2 = p2-O;

      /* calculate triangle edges */
      const avx3f e0 = v2-v0;
      const avx3f e1 = v0-v1;
      const avx3f e2 = v1-v2;

      /* calculate geometry normal and determinant */
      const avx3f Ng = cross(e1,e0);
      const avx3f Ng2 = Ng+Ng;
      const avxf det = dot(avx3f(Ng2),D);
      const avxf absDet = abs(det);
      const avxf sgnDet = signmsk(det);
      valid &= det != avxf(zero);
      if (unlikely(none(valid))) return;

      /* perform edge tests */
      const avxf U = dot(avx3f(cross(v2+v0,e0)),D) ^ sgnDet;
      valid &= U >= 0.0f;
      if (likely(none(valid))) return;
      const avxf V = dot(avx3f(cross(v0+v1,e1)),D) ^ sgnDet;
      valid &= V >= 0.0f;
      if (likely(none(valid))) return;
      const avxf W = dot(avx3f(cross(v1+v2,e2)),D) ^ sgnDet;
      valid &= W >= 0.0f;
      if (likely(none(valid))) return;
      
      /* perform depth test */
      const avxf T = dot(v0,avx3f(Ng2)) ^ sgnDet;
      valid &= (T >= absDet*ray.tnear) & (absDet*ray.tfar >= T);
      if (unlikely(none(valid))) return;

      /* update hit information for all rays that hit the triangle */
      const avxf rcpAbsDet = rcp(absDet);
      ray.u   = select(valid,U * rcpAbsDet,ray.u );
      ray.v   = select(valid,V * rcpAbsDet,ray.v );
      ray.tfar   = select(valid,T * rcpAbsDet,ray.tfar );
      ray.id0 = select(valid,id0,ray.id0);
      ray.id1 = select(valid,id1,ray.id1);
      ray.Ng.x = select(valid,Ng2.x,ray.Ng.x);
      ray.Ng.y = select(valid,Ng2.y,ray.Ng.y);
      ray.Ng.z = select(valid,Ng2.z,ray.Ng.z);
    }

    /*! Test for all rays if they are occluded by the triangle. */
    static __forceinline avxb occluded(const avxb& valid_i, const Ray8& ray, const Vector3f& p0, const Vector3f& p1, const Vector3f& p2)
    {
      /* calculate vertices relative to ray origin */
      avxb valid = valid_i;
      const avx3f O = ray.org;
      const avx3f D = ray.dir;
      const avx3f v0 = avx3f(p0)-O;
      const avx3f v1 = avx3f(p1)-O;
      const avx3f v2 = avx3f(p2)-O;

      /* calculate triangle edges */
      const avx3f e0 = v2-v0;
      const avx3f e1 = v0-v1;
      const avx3f e2 = v1-v2;

      /* calculate geometry normal and determinant */
      const avx3f Ng = cross(e1,e0);
      const avx3f Ng2 = Ng+Ng;
      const avxf det = dot(avx3f(Ng2),D);
      const avxf absDet = abs(det);
      const avxf sgnDet = signmsk(det);
      valid &= det != avxf(zero);
      if (unlikely(none(valid))) return valid;

      /* perform edge tests */
      const avxf U = dot(avx3f(cross(v2+v0,e0)),D) ^ sgnDet;
      valid &= U >= 0.0f;
      if (likely(none(valid))) return valid;
      const avxf V = dot(avx3f(cross(v0+v1,e1)),D) ^ sgnDet;
      valid &= V >= 0.0f;
      if (likely(none(valid))) return valid;
      const avxf W = dot(avx3f(cross(v1+v2,e2)),D) ^ sgnDet;
      valid &= W >= 0.0f;
      if (likely(none(valid))) return valid;
      
      /* perform depth test */
      const avxf T = dot(v0,avx3f(Ng2)) ^ sgnDet;
      valid &= (T >= absDet*ray.tnear) & (absDet*ray.tfar >= T);
      return valid;
    }

    /*! Test for all rays if they are occluded by the triangle. */
    static __forceinline avxb occluded(const avxb& valid_i, const Ray8& ray, const avx3f& p0, const avx3f& p1, const avx3f& p2)
    {
      /* calculate vertices relative to ray origin */
      avxb valid = valid_i;
      const avx3f O = ray.org;
      const avx3f D = ray.dir;
      const avx3f v0 = p0-O;
      const avx3f v1 = p1-O;
      const avx3f v2 = p2-O;

      /* calculate triangle edges */
      const avx3f e0 = v2-v0;
      const avx3f e1 = v0-v1;
      const avx3f e2 = v1-v2;

      /* calculate geometry normal and determinant */
      const avx3f Ng = cross(e1,e0);
      const avx3f Ng2 = Ng+Ng;
      const avxf det = dot(avx3f(Ng2),D);
      const avxf absDet = abs(det);
      const avxf sgnDet = signmsk(det);
      valid &= det != avxf(zero);
      if (unlikely(none(valid))) return valid;

      /* perform edge tests */
      const avxf U = dot(avx3f(cross(v2+v0,e0)),D) ^ sgnDet;
      valid &= U >= 0.0f;
      if (likely(none(valid))) return valid;
      const avxf V = dot(avx3f(cross(v0+v1,e1)),D) ^ sgnDet;
      valid &= V >= 0.0f;
      if (likely(none(valid))) return valid;
      const avxf W = dot(avx3f(cross(v1+v2,e2)),D) ^ sgnDet;
      valid &= W >= 0.0f;
      if (likely(none(valid))) return valid;
      
      /* perform depth test */
      const avxf T = dot(v0,avx3f(Ng2)) ^ sgnDet;
      valid &= (T >= absDet*ray.tnear) & (absDet*ray.tfar >= T);
      return valid;
    }
  };
}

#endif


