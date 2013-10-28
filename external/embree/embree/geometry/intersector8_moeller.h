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

#ifndef __EMBREE_INTERSECTOR8_MOELLER_H__
#define __EMBREE_INTERSECTOR8_MOELLER_H__

#include "intersector1_moeller.h"
#include "../common/ray8.h"

namespace embree
{
  /*! Implements a modified version of the Moeller
   *  Trumbore intersector for packets of 8 rays. */
  struct Intersector8MoellerTrumbore
  {
    typedef Intersector1MoellerTrumbore Intersector1;

    /*! Name of intersector */
    static const char* name() { return "moeller"; }

    /*! Intersect a ray with the 8 triangles and updates the hit. */
    static __forceinline void intersect(const avxb& valid_i, Ray8& ray, 
                                        const Vector3f& p0, const Vector3f& p1, const Vector3f& p2, const int& id0, const int& id1)
    {
      /* calculate edges and geometry normal */
      avxb valid = valid_i;
      const Vector3f e1 = p0-p1;
      const Vector3f e2 = p2-p0;
      const Vector3f Ng = cross(e1,e2);

      /* calculate determinant */
      const avx3f C = avx3f(p0) - ray.org;
      const avx3f R = cross(ray.dir,C);
      const avxf det = dot(avx3f(Ng),ray.dir);
      const avxf absDet = abs(det);
      const avxf sgnDet = signmsk(det);
      valid &= det != avxf(zero);
      if (unlikely(none(valid))) return;
      
      /* test against edge p2 p0 */
      const avxf U = dot(R,avx3f(e2)) ^ sgnDet;
      valid &= U >= 0.0f;
      if (likely(none(valid))) return;

      /* test against edge p0 p1 */
      const avxf V = dot(R,avx3f(e1)) ^ sgnDet;
      valid &= V >= 0.0f;
      if (likely(none(valid))) return;
      
      /* test against edge p1 p2 */
      const avxf W = absDet-U-V;
      valid &= W >= 0.0f;
      if (likely(none(valid))) return;

      /* perform depth test */
      const avxf T = dot(avx3f(Ng),C) ^ sgnDet;
      valid &= (T >= absDet*ray.tnear) & (absDet*ray.tfar >= T);
      if (unlikely(none(valid))) return;
      
      /* update hit information for all rays that hit the triangle */
      const avxf rcpAbsDet = rcp(absDet);
      ray.u   = select(valid,U * rcpAbsDet,ray.u );
      ray.v   = select(valid,V * rcpAbsDet,ray.v );
      ray.tfar   = select(valid,T * rcpAbsDet,ray.tfar );
      ray.id0 = select(valid,id0,ray.id0);
      ray.id1 = select(valid,id1,ray.id1);
      ray.Ng.x = select(valid,Ng.x,ray.Ng.x);
      ray.Ng.y = select(valid,Ng.y,ray.Ng.y);
      ray.Ng.z = select(valid,Ng.z,ray.Ng.z);
    }

    /*! Intersect a ray with the 8 triangles and updates the hit. */
    static __forceinline void intersect(const avxb& valid_i, Ray8& ray,
                                        const avx3f& p0, const avx3f& p1, const avx3f& p2, const avxi& id0, const avxi& id1)
    {
      /* calculate edges and geometry normal */
      avxb valid = valid_i;
      const avx3f e1 = p0-p1;
      const avx3f e2 = p2-p0;
      const avx3f Ng = cross(e1,e2);

      /* calculate determinant */
      const avx3f C = avx3f(p0) - ray.org;
      const avx3f R = cross(ray.dir,C);
      const avxf det = dot(avx3f(Ng),ray.dir);
      const avxf absDet = abs(det);
      const avxf sgnDet = signmsk(det);
      valid &= det != avxf(zero);
      if (unlikely(none(valid))) return;
      
      /* test against edge p2 p0 */
      const avxf U = dot(R,avx3f(e2)) ^ sgnDet;
      valid &= U >= 0.0f;
      if (likely(none(valid))) return;

      /* test against edge p0 p1 */
      const avxf V = dot(R,avx3f(e1)) ^ sgnDet;
      valid &= V >= 0.0f;
      if (likely(none(valid))) return;
      
      /* test against edge p1 p2 */
      const avxf W = absDet-U-V;
      valid &= W >= 0.0f;
      if (likely(none(valid))) return;

      /* perform depth test */
      const avxf T = dot(avx3f(Ng),C) ^ sgnDet;
      valid &= (T >= absDet*ray.tnear) & (absDet*ray.tfar >= T);
      if (unlikely(none(valid))) return;
      
      /* update hit information for all rays that hit the triangle */
      const avxf rcpAbsDet = rcp(absDet);
      ray.u   = select(valid,U * rcpAbsDet,ray.u );
      ray.v   = select(valid,V * rcpAbsDet,ray.v );
      ray.tfar = select(valid,T * rcpAbsDet,ray.tfar );
      ray.id0 = select(valid,id0,ray.id0);
      ray.id1 = select(valid,id1,ray.id1);
      ray.Ng.x = select(valid,Ng.x,ray.Ng.x);
      ray.Ng.y = select(valid,Ng.y,ray.Ng.y);
      ray.Ng.z = select(valid,Ng.z,ray.Ng.z);
    }

    /*! Test for all rays if they are occluded by the triangle. */
    static __forceinline avxb occluded(const avxb& valid_i, const Ray8& ray, const Vector3f& p0, const Vector3f& p1, const Vector3f& p2)
    {
      /* calculate edges and geometry normal */
      avxb valid = valid_i;
      const Vector3f e1 = p0-p1;
      const Vector3f e2 = p2-p0;
      const Vector3f Ng = cross(e1,e2);

      /* calculate determinant */
      const avx3f C = avx3f(p0) - ray.org;
      const avx3f R = cross(ray.dir,C);
      const avxf det = dot(avx3f(Ng),ray.dir);
      const avxf absDet = abs(det);
      const avxf sgnDet = signmsk(det);
      valid &= det != avxf(zero);
      if (unlikely(none(valid))) return valid;
      
      /* test against edge p2 p0 */
      const avxf U = dot(R,avx3f(e2)) ^ sgnDet;
      valid &= U >= 0.0f;
      if (likely(none(valid))) return valid;

      /* test against edge p0 p1 */
      const avxf V = dot(R,avx3f(e1)) ^ sgnDet;
      valid &= V >= 0.0f;
      if (likely(none(valid))) return valid;
      
      /* test against edge p1 p2 */
      const avxf W = absDet-U-V;
      valid &= W >= 0.0f;
      if (likely(none(valid))) return valid;

      /* perform depth test */
      const avxf T = dot(avx3f(Ng),C) ^ sgnDet;
      valid &= (T >= absDet*ray.tnear) & (absDet*ray.tfar >= T);
      return valid;
    }

    /*! Test for all rays if they are occluded by the triangle. */
    static __forceinline avxb occluded(const avxb& valid_i, const Ray8& ray, const avx3f& p0, const avx3f& p1, const avx3f& p2)
    {
      /* calculate edges and geometry normal */
      avxb valid = valid_i;
      const avx3f e1 = p0-p1;
      const avx3f e2 = p2-p0;
      const avx3f Ng = cross(e1,e2);

      /* calculate determinant */
      const avx3f C = avx3f(p0) - ray.org;
      const avx3f R = cross(ray.dir,C);
      const avxf det = dot(avx3f(Ng),ray.dir);
      const avxf absDet = abs(det);
      const avxf sgnDet = signmsk(det);
      valid &= det != avxf(zero);
      if (unlikely(none(valid))) return valid;
      
      /* test against edge p2 p0 */
      const avxf U = dot(R,avx3f(e2)) ^ sgnDet;
      valid &= U >= 0.0f;
      if (likely(none(valid))) return valid;

      /* test against edge p0 p1 */
      const avxf V = dot(R,avx3f(e1)) ^ sgnDet;
      valid &= V >= 0.0f;
      if (likely(none(valid))) return valid;
      
      /* test against edge p1 p2 */
      const avxf W = absDet-U-V;
      valid &= W >= 0.0f;
      if (likely(none(valid))) return valid;

      /* perform depth test */
      const avxf T = dot(avx3f(Ng),C) ^ sgnDet;
      valid &= (T >= absDet*ray.tnear) & (absDet*ray.tfar >= T);
      return valid;
    }
  };
}

#endif

