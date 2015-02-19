// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
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

#pragma once

#include "triangle1v.h"
#include "common/ray.h"
#include "geometry/filter.h"

namespace embree
{
  namespace isa
  {
    /*! Modified Pluecker ray/triangle intersector. The test first shifts the ray
     *  origin into the origin of the coordinate system and then uses
     *  Pluecker coordinates for the intersection. Due to the shift, the
     *  Pluecker coordinate calculation simplifies. The edge equations
     *  are watertight along the edge for neighboring triangles. */
    template<bool list>
      struct Triangle1vIntersector1Pluecker
      {
        typedef Triangle1v Primitive;
        
        struct Precalculations {
          __forceinline Precalculations (const Ray& ray) {}
        };
        
        /*! Intersect a ray with the triangle and updates the hit. */
        static __forceinline void intersect(const Precalculations& pre, Ray& ray, const Primitive& tri, Scene* scene)
        {
          /* load triangle */
          STAT3(normal.trav_prims,1,1,1);
          const Vec3fa tri_v0 = tri.v0;
          const Vec3fa tri_v1 = tri.v1;
          const Vec3fa tri_v2 = tri.v2;
          
          /* calculate vertices relative to ray origin */
          const Vec3fa O = ray.org;
          const Vec3fa D = ray.dir;
          const Vec3fa v0 = tri_v0-O;
          const Vec3fa v1 = tri_v1-O;
          const Vec3fa v2 = tri_v2-O;
          
          /* calculate triangle edges */
          const Vec3fa e0 = v2-v0;
          const Vec3fa e1 = v0-v1;
          const Vec3fa e2 = v1-v2;
          
          /* calculate geometry normal and denominator */
          const Vec3fa Ng1 = cross(e1,e0);
          const Vec3fa Ng = Ng1+Ng1;
          const float den = dot(Ng,D);
          const float absDen = abs(den);
          const float sgnDen = signmsk(den);
          
          /* perform edge tests */
          const float U = xorf(dot(cross(v2+v0,e0),D),sgnDen);
          if (unlikely(U < 0.0f)) return;
          const float V = xorf(dot(cross(v0+v1,e1),D),sgnDen);
          if (unlikely(V < 0.0f)) return;
          const float W = xorf(dot(cross(v1+v2,e2),D),sgnDen);
          if (unlikely(W < 0.0f)) return;
          
          /* perform depth test */
          const float T = xorf(dot(v0,Ng),sgnDen);
          if (unlikely(absDen*float(ray.tfar) < T)) return;
          if (unlikely(T < absDen*float(ray.tnear))) return;
          
          /* perform backface culling */
#if defined(RTCORE_BACKFACE_CULLING)
          if (unlikely(den <= 0.0f)) return;
#else
          if (unlikely(den == 0.0f)) return;
#endif
          
          /* ray masking test */
#if defined(RTCORE_RAY_MASK)
          if (unlikely((tri.mask() & ray.mask) == 0)) return;
#endif
          
          /* calculate hit information */
          const float rcpAbsDen = rcp(absDen);
          const float u = U * rcpAbsDen;
          const float v = V * rcpAbsDen;
          const float t = T * rcpAbsDen;
          const int geomID = tri.geomID<list>();
          const int primID = tri.primID<list>();
          
          /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER)
          Geometry* geometry = scene->get(geomID);
          if (unlikely(geometry->hasIntersectionFilter1())) {
            runIntersectionFilter1(geometry,ray,u,v,t,Ng,geomID,primID);
            return;
          }
#endif
          
          /* update hit information */
          ray.u = u;
          ray.v = v;
          ray.tfar = t;
          ray.Ng  = Ng;
          ray.geomID = geomID;
          ray.primID = primID;
        }
        
        /*! Test if the ray is occluded by one of the triangles. */
        static __forceinline bool occluded(const Precalculations& pre, Ray& ray, const Primitive& tri, Scene* scene)
        {
          /* load triangle */
          STAT3(shadow.trav_prims,1,1,1);
          const Vec3fa tri_v0 = tri.v0;
          const Vec3fa tri_v1 = tri.v1;
          const Vec3fa tri_v2 = tri.v2;
          
          /* calculate vertices relative to ray origin */
          const Vec3fa O = ray.org;
          const Vec3fa D = ray.dir;
          const Vec3fa v0 = tri_v0-O;
          const Vec3fa v1 = tri_v1-O;
          const Vec3fa v2 = tri_v2-O;
          
          /* calculate triangle edges */
          const Vec3fa e0 = v2-v0;
          const Vec3fa e1 = v0-v1;
          const Vec3fa e2 = v1-v2;
          
          /* calculate geometry normal and denominator */
          const Vec3fa Ng1 = cross(e1,e0);
          const Vec3fa Ng = Ng1+Ng1;
          const float den = dot(Ng,D);
          const float absDen = abs(den);
          const float sgnDen = signmsk(den);
          
          /* perform edge tests */
          const float U = xorf(dot(cross(v2+v0,e0),D),sgnDen);
          if (unlikely(U < 0.0f)) return false;
          const float V = xorf(dot(cross(v0+v1,e1),D),sgnDen);
          if (unlikely(V < 0.0f)) return false;
          const float W = xorf(dot(cross(v1+v2,e2),D),sgnDen);
          if (unlikely(W < 0.0f)) return false;
          
          /* perform depth test */
          const float T = xorf(dot(v0,Ng),sgnDen);
          if (unlikely(absDen*float(ray.tfar) < T)) return false;
          if (unlikely(T < absDen*float(ray.tnear))) return false;
          
          /* perform backface culling */
#if defined(RTCORE_BACKFACE_CULLING)
          if (unlikely(den <= 0.0f)) return false;
#else
          if (unlikely(den == 0.0f)) return false;
#endif
          
          /* ray masking test */
#if defined(RTCORE_RAY_MASK)
          if (unlikely((tri.mask() & ray.mask) == 0)) return false;
#endif
          
          /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER)
          const int geomID = tri.geomID<list>();
          Geometry* geometry = scene->get(geomID);
          if (unlikely(geometry->hasOcclusionFilter1()))
          {
            /* calculate hit information */
            const float rcpAbsDen = rcp(absDen);
            const float u = U*rcpAbsDen;
            const float v = V*rcpAbsDen;
            const float t = T*rcpAbsDen;
            const int primID = tri.primID<list>();
            return runOcclusionFilter1(geometry,ray,u,v,t,Ng,geomID,primID);
          }
#endif
          
          return true;
        }
      };
  }
}
