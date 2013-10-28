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

#ifndef __EMBREE_INTERSECTOR1_PLUECKER_H__
#define __EMBREE_INTERSECTOR1_PLUECKER_H__

#include "../common/ray.h"

namespace embree
{
  /*! Modified Pluecker ray/triangle intersector. The test first shifts the ray
   *  origin into the origin of the coordinate system and then uses
   *  Pluecker coordinates for the intersection. Due to the shift, the
   *  Pluecker coordinate calculation simplifies. The edge equations
   *  are watertight along the edge for neighboring triangles. */
  struct Intersector1Pluecker
  {
    /*! Name of intersector */
    static const char* name() { return "pluecker"; }

    /*! Intersect a ray with the 4 triangles and updates the hit. */
    static __forceinline void intersect(Ray& ray, const Vector3f& p0, const Vector3f& p1, const Vector3f& p2, const int& id0, const int& id1)
    {
      /* calculate vertices relative to ray origin */
      const Vector3f O = ray.org;
      const Vector3f D = ray.dir;
      const Vector3f v0 = p0-O;
      const Vector3f v1 = p1-O;
      const Vector3f v2 = p2-O;

      /* calculate triangle edges */
      const Vector3f e0 = v2-v0;
      const Vector3f e1 = v0-v1;
      const Vector3f e2 = v1-v2;

      /* calculate geometry normal and determinant */
      const Vector3f Ng = cross(e1,e0);
      const Vector3f Ng2 = Ng+Ng;
      const float det = dot(Ng2,D);
      const float absDet = abs(det);
      const float sgnDet = signmsk(det);

      /* perform edge tests */
      const float U = xorf(dot(cross(v2+v0,e0),D),sgnDet);
      if (unlikely(U < 0.0f)) return;
      const float V = xorf(dot(cross(v0+v1,e1),D),sgnDet);
      if (unlikely(V < 0.0f)) return;
      const float W = xorf(dot(cross(v1+v2,e2),D),sgnDet);
      if (unlikely(W < 0.0f)) return;
      
      /* perform depth test */
      const float T = xorf(dot(v0,Ng2),sgnDet);
      if (unlikely(absDet*float(ray.tfar) < T)) return;
      if (unlikely(T < absDet*float(ray.tnear))) return;
      if (unlikely(det == float(zero))) return;

      /* update hit information */
      const float rcpAbsDet = rcp(absDet);
      ray.u   = U * rcpAbsDet;
      ray.v   = V * rcpAbsDet;
      ray.tfar = T * rcpAbsDet;
      ray.Ng  = Ng2;
      ray.id0 = id0;
      ray.id1 = id1;
    }

    /*! Test if the ray is occluded by one of the triangles. */
    static __forceinline bool occluded(const Ray& ray, const Vector3f& p0, const Vector3f& p1, const Vector3f& p2)
    {
      /* calculate vertices relative to ray origin */
      const Vector3f O = ray.org;
      const Vector3f D = ray.dir;
      const Vector3f v0 = p0-O;
      const Vector3f v1 = p1-O;
      const Vector3f v2 = p2-O;

      /* calculate triangle edges */
      const Vector3f e0 = v2-v0;
      const Vector3f e1 = v0-v1;
      const Vector3f e2 = v1-v2;

      /* calculate geometry normal and determinant */
      const Vector3f Ng = cross(e1,e0);
      const Vector3f Ng2 = Ng+Ng;
      const float det = dot(Ng2,D);
      const float absDet = abs(det);
      const float sgnDet = signmsk(det);

      /* perform edge tests */
      const float U = xorf(dot(cross(v2+v0,e0),D),sgnDet);
      if (unlikely(U < 0.0f)) return false;
      const float V = xorf(dot(cross(v0+v1,e1),D),sgnDet);
      if (unlikely(V < 0.0f)) return false;
      const float W = xorf(dot(cross(v1+v2,e2),D),sgnDet);
      if (unlikely(W < 0.0f)) return false;
      
      /* perform depth test */
      const float T = xorf(dot(v0,Ng2),sgnDet);
      if (unlikely(absDet*float(ray.tfar) < T)) return false;
      if (unlikely(T < absDet*float(ray.tnear))) return false;
      if (unlikely(det == float(zero))) return false;
      return true;
    }

#if defined (__SSE__)

    /*! Intersect a ray with the 4 triangles and updates the hit. */
    static __forceinline void intersect(Ray& ray, const sse3f& p0, const sse3f& p1, const sse3f& p2, const ssei& id0, const ssei& id1)
    {
      /* calculate vertices relative to ray origin */
      const sse3f O = sse3f(ray.org);
      const sse3f D = sse3f(ray.dir);
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
      const ssef det = dot(Ng2,D);
      const ssef absDet = abs(det);
      const ssef sgnDet = signmsk(det);

      /* perform edge tests */
      const ssef U = dot(cross(v2+v0,e0),D) ^ sgnDet;
      const ssef V = dot(cross(v0+v1,e1),D) ^ sgnDet;
      const ssef W = dot(cross(v1+v2,e2),D) ^ sgnDet;
      sseb valid = (det != ssef(zero)) & (U >= 0.0f) & (V >= 0.0f) & (W >= 0.0f);
      if (unlikely(none(valid))) return;

      /* perform depth test */
      const ssef T = dot(v0,Ng2) ^ sgnDet;
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
      ray.Ng.x = Ng2.x[i];
      ray.Ng.y = Ng2.y[i];
      ray.Ng.z = Ng2.z[i];
      ray.id0 = id0[i] & 0x7FFFFFFF;
      ray.id1 = id1[i];
    }

    /*! Test if the ray is occluded by one of the triangles. */
    static __forceinline bool occluded(const Ray& ray, const sse3f& p0, const sse3f& p1, const sse3f& p2)
    {
      STAT3(shadow.trav_tris,1,1,1);

      /* calculate vertices relative to ray origin */
      const sse3f O = sse3f(ray.org);
      const sse3f D = sse3f(ray.dir);
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
      const ssef det = dot(Ng2,D);
      const ssef absDet = abs(det);
      const ssef sgnDet = signmsk(det);

      /* perform edge tests */
      const ssef U = dot(cross(v2+v0,e0),D) ^ sgnDet;
      const ssef V = dot(cross(v0+v1,e1),D) ^ sgnDet;
      const ssef W = dot(cross(v1+v2,e2),D) ^ sgnDet;
      sseb valid = (det != ssef(zero)) & (U >= 0.0f) & (V >= 0.0f) & (W >= 0.0f);
      if (unlikely(none(valid))) return false;
      
      /* perform depth test */
      const ssef T = dot(v0,Ng2) ^ sgnDet;
      valid &= (T >= absDet*ssef(ray.tnear)) & (absDet*ssef(ray.tfar) >= T);
      if (unlikely(none(valid))) return false;

      return true;
    }
#endif
  };
}

#endif
