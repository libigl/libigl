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

#ifndef __EMBREE_ACCEL_TRIANGLE1I_INTERSECTOR_PLUECKER_H__
#define __EMBREE_ACCEL_TRIANGLE1I_INTERSECTOR_PLUECKER_H__

#include "triangle1.h"
#include "../common/ray.h"
#include "../common/hit.h"

namespace embree
{
  /*! Modified Pluecker ray/triangle intersector for a single ray with
   *  individual indexed triangles. The test first shifts the ray
   *  origin into the origin of the coordinate system and then uses
   *  Pluecker coordinates for the intersection. Due to the shift, the
   *  Pluecker coordinate calculation simplifies. The edge equations
   *  are watertight along the edge for neighboring triangles. */
  struct Triangle1iIntersectorPluecker
  {
    typedef Triangle1i Triangle;

    /*! Intersect a ray with the 4 triangles and updates the hit. */
    static __forceinline void intersect(const Ray& ray, Hit& hit, const Triangle1i& tri, const Vec3fa* vertices)
    {
      STAT3(normal.trav_tris,1,1,1);

      /* calculate edges, geometry normal, and determinant */
      const Vec3f O = ray.org;
      const Vec3f D = ray.dir;
      const Vec3f v0 = vertices[tri.v0]-O;
      const Vec3f v1 = vertices[tri.v1]-O;
      const Vec3f v2 = vertices[tri.v2]-O;
      const Vec3f e0 = v2-v0;
      const Vec3f e1 = v0-v1;
      const Vec3f e2 = v1-v2;
      const Vec3f Ng = 2.0f*cross(e1,e0);
      const float det = dot(Ng,D);
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
      const float T = xorf(dot(v0,Ng),sgnDet);
      if (unlikely(absDet*float(hit.t) < T)) return;
      if (unlikely(T < absDet*float(ray.near))) return;
      if (unlikely(det == float(zero))) return;

      /* update hit information */
      const float rcpAbsDet = rcp(absDet);
      hit.u   = U * rcpAbsDet;
      hit.v   = V * rcpAbsDet;
      hit.t   = T * rcpAbsDet;
      hit.id0 = tri.id0;
      hit.id1 = tri.id1;
    }

    /*! Test if the ray is occluded by one of the triangles. */
    static __forceinline bool occluded(const Ray& ray, const Triangle1i& tri, const Vec3fa* vertices = NULL)
    {
      STAT3(shadow.trav_tris,1,1,1);

      /* calculate edges, geometry normal, and determinant */
      const Vec3f O = ray.org;
      const Vec3f D = ray.dir;
      const Vec3f v0 = vertices[tri.v0]-O;
      const Vec3f v1 = vertices[tri.v1]-O;
      const Vec3f v2 = vertices[tri.v2]-O;
      const Vec3f e0 = v2-v0;
      const Vec3f e1 = v0-v1;
      const Vec3f e2 = v1-v2;
      const Vec3f Ng = 2.0f*cross(e1,e0);
      const float det = dot(Ng,D);
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
      const float T = xorf(dot(v0,Ng),sgnDet);
      if (unlikely(absDet*float(ray.far) < T)) return false;
      if (unlikely(T < absDet*float(ray.near))) return false;
      if (unlikely(det == float(zero))) return false;
      return true;
    }
  };
}

#endif


