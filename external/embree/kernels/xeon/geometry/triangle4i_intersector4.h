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

#ifndef __EMBREE_TRIANGLE4I_INTERSECTOR4_PLUECKER_H__
#define __EMBREE_TRIANGLE4I_INTERSECTOR4_PLUECKER_H__

#include "triangle4i.h"
#include "common/ray4.h"

namespace embree
{
  /*! Intersector4 for triangle4i */
  struct Triangle4iIntersector4Pluecker
  {
    typedef Triangle4i Primitive;

    static __forceinline void intersect(const sseb& valid_i, Ray4& ray, const Triangle4i& tri, const void* geom)
    {
      for (size_t i=0; i<tri.size(); i++)
      {
        STAT3(normal.trav_prims,1,popcnt(valid_i),4);

        /* load vertices */
        const Vec3f& p0 = *tri.v0[i];
        const Vec3f& p1 = *(Vec3f*)((int*)&p0 + tri.v1[i]);
        const Vec3f& p2 = *(Vec3f*)((int*)&p0 + tri.v2[i]);

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
        
        /* calculate geometry normal and denominator */
        const sse3f Ng1 = cross(e1,e0);
        const sse3f Ng = Ng1+Ng1;
        const ssef den = dot(sse3f(Ng),D);
        const ssef absDen = abs(den);
        const ssef sgnDen = signmsk(den);
        
        /* perform edge tests */
        const ssef U = dot(sse3f(cross(v2+v0,e0)),D) ^ sgnDen;
        valid &= U >= 0.0f;
        if (likely(none(valid))) continue;
        const ssef V = dot(sse3f(cross(v0+v1,e1)),D) ^ sgnDen;
        valid &= V >= 0.0f;
        if (likely(none(valid))) continue;
        const ssef W = dot(sse3f(cross(v1+v2,e2)),D) ^ sgnDen;
        valid &= W >= 0.0f;
        if (likely(none(valid))) continue;
        
        /* perform depth test */
        const ssef T = dot(v0,sse3f(Ng)) ^ sgnDen;
        valid &= (T >= absDen*ray.tnear) & (absDen*ray.tfar >= T);
        if (unlikely(none(valid))) continue;

        /* perform backface culling */
#if defined(__BACKFACE_CULLING__)
        valid &= den > ssef(zero);
        if (unlikely(none(valid))) continue;
#else
        valid &= den != ssef(zero);
        if (unlikely(none(valid))) continue;
#endif

        /* ray masking test */
#if defined(__USE_RAY_MASK__)
        int mask = ((Scene*)geom)->getTriangleMesh(tri.geomID[i])->mask;
        valid &= (mask & ray.mask) != 0;
        if (unlikely(none(valid))) continue;
#endif
        
        /* calculate hit information */
        const ssef u = U / absDen;
        const ssef v = V / absDen;
        const ssef t = T / absDen;
        const int geomID = tri.geomID[i];
        const int primID = tri.primID[i];

        /* intersection filter test */
#if defined(__INTERSECTION_FILTER__)
        Geometry* geometry = ((Scene*)geom)->get(geomID);
        if (unlikely(geometry->hasIntersectionFilter4())) {
          runIntersectionFilter4(valid,geometry,ray,u,v,t,Ng,geomID,primID);
          continue;
        }
#endif

        /* update hit information */
        store4f(valid,&ray.u,u);
        store4f(valid,&ray.v,v);
        store4f(valid,&ray.tfar,t);
        store4i(valid,&ray.geomID,geomID);
        store4i(valid,&ray.primID,primID);
        store4f(valid,&ray.Ng.x,Ng.x);
        store4f(valid,&ray.Ng.y,Ng.y);
        store4f(valid,&ray.Ng.z,Ng.z);
      }
    }

    static __forceinline void intersect(const sseb& valid, Ray4& ray, const Triangle4i* tri, size_t num, const void* geom)
    {
      for (size_t i=0; i<num; i++)
        intersect(valid,ray,tri[i],geom);
    }
    
    static __forceinline sseb occluded(const sseb& valid_i, Ray4& ray, const Triangle4i& tri, const void* geom)
    {
      sseb valid0 = valid_i;

      for (size_t i=0; i<tri.size(); i++)
      {
        STAT3(shadow.trav_prims,1,popcnt(valid_i),4);

        /* load vertices */
        const Vec3f& p0 = *tri.v0[i];
        const Vec3f& p1 = *(Vec3f*)((int*)&p0 + tri.v1[i]);
        const Vec3f& p2 = *(Vec3f*)((int*)&p0 + tri.v2[i]);

        /* calculate vertices relative to ray origin */
        sseb valid = valid0;
        const sse3f O = ray.org;
        const sse3f D = ray.dir;
        const sse3f v0 = sse3f(p0)-O;
        const sse3f v1 = sse3f(p1)-O;
        const sse3f v2 = sse3f(p2)-O;

        /* calculate triangle edges */
        const sse3f e0 = v2-v0;
        const sse3f e1 = v0-v1;
        const sse3f e2 = v1-v2;
        
        /* calculate geometry normal and denominator */
        const sse3f Ng1 = cross(e1,e0);
        const sse3f Ng = Ng1+Ng1;
        const ssef den = dot(sse3f(Ng),D);
        const ssef absDen = abs(den);
        const ssef sgnDen = signmsk(den);
        
        /* perform edge tests */
        const ssef U = dot(sse3f(cross(v2+v0,e0)),D) ^ sgnDen;
        valid &= U >= 0.0f;
        if (likely(none(valid))) continue;
        const ssef V = dot(sse3f(cross(v0+v1,e1)),D) ^ sgnDen;
        valid &= V >= 0.0f;
        if (likely(none(valid))) continue;
        const ssef W = dot(sse3f(cross(v1+v2,e2)),D) ^ sgnDen;
        valid &= W >= 0.0f;
        if (likely(none(valid))) continue;
        
        /* perform depth test */
        const ssef T = dot(v0,sse3f(Ng)) ^ sgnDen;
        valid &= (T >= absDen*ray.tnear) & (absDen*ray.tfar >= T);

        /* perform backface culling */
#if defined(__BACKFACE_CULLING__)
        valid &= den > ssef(zero);
        if (unlikely(none(valid))) continue;
#else
        valid &= den != ssef(zero);
        if (unlikely(none(valid))) continue;
#endif

        /* ray masking test */
#if defined(__USE_RAY_MASK__)
        int mask = ((Scene*)geom)->getTriangleMesh(tri.geomID[i])->mask;
        valid &= (mask & ray.mask) != 0;
        if (unlikely(none(valid))) continue;
#endif

        /* intersection filter test */
#if defined(__INTERSECTION_FILTER__)
        const int geomID = tri.geomID[i];
        Geometry* geometry = ((Scene*)geom)->get(geomID);
        if (unlikely(geometry->hasOcclusionFilter4()))
        {
          /* calculate hit information */
          const ssef u = U / absDen;
          const ssef v = V / absDen;
          const ssef t = T / absDen;
          const int primID = tri.primID[i];
          valid = runOcclusionFilter4(valid,geometry,ray,u,v,t,Ng,geomID,primID);
        }
#endif

        /* update occlusion */
        valid0 &= !valid;
        if (none(valid0)) break;
      }
      return !valid0;
    }

    static __forceinline sseb occluded(const sseb& valid, Ray4& ray, const Triangle4i* tri, size_t num, const void* geom)
    {
      sseb valid0 = valid;
      for (size_t i=0; i<num; i++) {
        valid0 &= !occluded(valid0,ray,tri[i],geom);
        if (none(valid0)) break;
      }
      return !valid0;
    }
  };
}

#endif


