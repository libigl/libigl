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

#ifndef __EMBREE_INTERSECTOR1_MOELLER_H__
#define __EMBREE_INTERSECTOR1_MOELLER_H__

#include "../common/ray.h"

namespace embree
{
  /*! Implements a modified version of the Moeller Trumbore intersector. */
  struct Intersector1MoellerTrumbore
  {
    /*! Name of intersector */
    static const char* name() { return "moeller"; }

    /*! Intersect a single ray with a single triangle. */
    static __forceinline void intersect(Ray& ray, const Vector3f& p0, const Vector3f& p1, const Vector3f& p2, const int& id0, const int& id1)
    {
      /* calculate edges and geometry normal */
      const Vector3f e1 = p0-p1;
      const Vector3f e2 = p2-p0;
      const Vector3f Ng = cross(e1,e2);

      /* calculate determinant */
      const Vector3f C = p0 - ray.org;
      const Vector3f R = cross(ray.dir,C);
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
      if (unlikely(T < absDet*ray.tnear || absDet*ray.tfar < T)) return;

      /* update hit information */
      const float rcpAbsDet = rcp(absDet);
      ray.u   = U * rcpAbsDet;
      ray.v   = V * rcpAbsDet;
      ray.tfar = T * rcpAbsDet;
      ray.Ng  = Ng;
      ray.id0 = id0;
      ray.id1 = id1;
    }

    /*! Test if the ray is occluded by one of the triangles. */
    static __forceinline bool occluded(const Ray& ray, const Vector3f& p0, const Vector3f& p1, const Vector3f& p2)
    {
      /* calculate edges and geometry normal */
      const Vector3f e1 = p0-p1;
      const Vector3f e2 = p2-p0;
      const Vector3f Ng = cross(e1,e2);
      
      /* calculate determinant */
      const Vector3f C = p0 - ray.org;
      const Vector3f R = cross(ray.dir,C);
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
      if (unlikely(T < absDet*ray.tnear || absDet*ray.tfar < T)) return false;

      return true;
    }

    /*! Intersect a ray with the 4 triangles and updates the hit. */
    template<typename Ray>
      static __forceinline void intersect(const size_t k, Ray& ray, const Vector3f& p0, const Vector3f& p1, const Vector3f& p2, const int& id0, const int& id1)
    {
      /* calculate edges and geometry normal */
      const Vector3f org(ray.org.x[k],ray.org.y[k],ray.org.z[k]);
      const Vector3f dir(ray.dir.x[k],ray.dir.y[k],ray.dir.z[k]);
      const Vector3f e1 = p0-p1;
      const Vector3f e2 = p2-p0;
      const Vector3f Ng = cross(e1,e2);

      /* calculate determinant */
      const Vector3f C = p0 - org;
      const Vector3f R = cross(dir,C);
      const float det = dot(Ng,dir);
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
      if (unlikely(T < absDet*ray.tnear[k] || absDet*ray.tfar[k] < T)) return;

      /* update hit information */
      const float rcpAbsDet = rcp(absDet);
      ray.u[k]   = U * rcpAbsDet;
      ray.v[k]   = V * rcpAbsDet;
      ray.tfar[k] = T * rcpAbsDet;
      ray.Ng.x[k]  = Ng.x;
      ray.Ng.y[k]  = Ng.y;
      ray.Ng.z[k]  = Ng.z;
      ray.id0[k] = id0;
      ray.id1[k] = id1;
    }

    /*! Test if the ray is occluded by one of the triangles. */
    template<typename Ray>
      static __forceinline bool occluded(const size_t k, const Ray& ray, const Vector3f& p0, const Vector3f& p1, const Vector3f& p2)
    {
      /* calculate edges and geometry normal */
      const Vector3f org(ray.org.x[k],ray.org.y[k],ray.org.z[k]);
      const Vector3f dir(ray.dir.x[k],ray.dir.y[k],ray.dir.z[k]);
      const Vector3f e1 = p0-p1;
      const Vector3f e2 = p2-p0;
      const Vector3f Ng = cross(e1,e2);
      
      /* calculate determinant */
      const Vector3f C = p0 - org;
      const Vector3f R = cross(dir,C);
      const float det = dot(Ng,dir);
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
      if (unlikely(T < absDet*ray.tnear[k] || absDet*ray.tfar[k] < T)) return false;

      return true;
    }
    
#if defined (__SSE__)

    /*! Intersect a ray with the 4 triangles and updates the hit. */
    static __forceinline void intersect(Ray& ray, const sse3f& p0, const sse3f& p1, const sse3f& p2, const ssei& id0, const ssei& id1)
    {
      /* calculate edges and geometry normal */
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
      sseb valid = (det != ssef(zero)) & (U >= 0.0f) & (V >= 0.0f) & (W >= 0.0f);
      if (unlikely(none(valid))) return;
      
      /* perform depth test */
      const ssef T = dot(Ng,C) ^ sgnDet;
      valid &= (T >= absDet*ssef(ray.tnear)) & (absDet*ssef(ray.tfar) >= T);
      if (unlikely(none(valid))) return;

      /* update hit information */
      const ssef rcpAbsDet = rcp(absDet);
      const ssef u = U * rcpAbsDet;
      const ssef v = V * rcpAbsDet;
      const ssef t = T * rcpAbsDet;
      const size_t i = select_min(valid,t);
      ray.tfar = t[i];
      ray.u = u[i];
      ray.v = v[i];
      ray.Ng.x = Ng.x[i];
      ray.Ng.y = Ng.y[i];
      ray.Ng.z = Ng.z[i];
      ray.id0 = id0[i];
      ray.id1 = id1[i];
    }

    /*! Test if the ray is occluded by one of the triangles. */
    static __forceinline bool occluded(const Ray& ray, const sse3f& p0, const sse3f& p1, const sse3f& p2)
    {
      /* calculate edges and geometry normal */
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
      sseb valid = (det != ssef(zero)) & (U >= 0.0f) & (V >= 0.0f) & (W >= 0.0f);
      if (unlikely(none(valid))) return false;
      
      /* perform depth test */
      const ssef T = dot(Ng,C) ^ sgnDet;
      valid &= (T >= absDet*ssef(ray.tnear)) & (absDet*ssef(ray.tfar) >= T);
      if (unlikely(none(valid))) return false;

      return true;
    }

#endif
  };
}

#endif


