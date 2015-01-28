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

#ifndef __EMBREE_TRIANGLE4I_INTERSECTOR8_PLUECKER_H__
#define __EMBREE_TRIANGLE4I_INTERSECTOR8_PLUECKER_H__

#include "triangle4i.h"
#include "common/ray4.h"

namespace embree
{
  /*! Intersector8 for triangle4i */
  struct Triangle4iIntersector8Pluecker
  {
    typedef Triangle4i Primitive;
  
    static __forceinline void intersect(const avxb& valid_i, Ray8& ray, const Triangle4i& tri, const void* geom)
    {
      for (size_t i=0; i<tri.size(); i++)
      {
        STAT3(normal.trav_prims,1,popcnt(valid_i),8);

        /* load vertices */
        const Vec3f& p0 = *tri.v0[i];
        const Vec3f& p1 = *(Vec3f*)((int*)&p0 + tri.v1[i]);
        const Vec3f& p2 = *(Vec3f*)((int*)&p0 + tri.v2[i]);

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
        
        /* calculate geometry normal and denominator */
        const avx3f Ng1 = cross(e1,e0);
        const avx3f Ng = Ng1+Ng1;
        const avxf den = dot(avx3f(Ng),D);
        const avxf absDen = abs(den);
        const avxf sgnDen = signmsk(den);
        
        /* perform edge tests */
        const avxf U = dot(avx3f(cross(v2+v0,e0)),D) ^ sgnDen;
        valid &= U >= 0.0f;
        if (likely(none(valid))) continue;
        const avxf V = dot(avx3f(cross(v0+v1,e1)),D) ^ sgnDen;
        valid &= V >= 0.0f;
        if (likely(none(valid))) continue;
        const avxf W = dot(avx3f(cross(v1+v2,e2)),D) ^ sgnDen;
        valid &= W >= 0.0f;
        if (likely(none(valid))) continue;
        
        /* perform depth test */
        const avxf T = dot(v0,avx3f(Ng)) ^ sgnDen;
        valid &= (T >= absDen*ray.tnear) & (absDen*ray.tfar >= T);
        if (unlikely(none(valid))) continue;

        /* perform backface culling */
#if defined(__BACKFACE_CULLING__)
        valid &= den > avxf(zero);
        if (unlikely(none(valid))) continue;
#else
        valid &= den != avxf(zero);
        if (unlikely(none(valid))) continue;
#endif

        /* ray masking test */
#if defined(__USE_RAY_MASK__)
        int mask = ((Scene*)geom)->getTriangleMesh(tri.geomID[i])->mask;
        valid &= (mask & ray.mask) != 0;
        if (unlikely(none(valid))) continue;
#endif
        
        /* calculate hit information */
        const avxf u = U / absDen;
        const avxf v = V / absDen;
        const avxf t = T / absDen;
        const int geomID = tri.geomID[i];
        const int primID = tri.primID[i];

        /* intersection filter test */
#if defined(__INTERSECTION_FILTER__)
        Geometry* geometry = ((Scene*)geom)->get(geomID);
        if (unlikely(geometry->hasIntersectionFilter8())) {
          runIntersectionFilter8(valid,geometry,ray,u,v,t,Ng,geomID,primID);
          continue;
        }
#endif

        /* update hit information */
        store8f(valid,&ray.u,u);
        store8f(valid,&ray.v,v);
        store8f(valid,&ray.tfar,t);
        store8i(valid,&ray.geomID,geomID);
        store8i(valid,&ray.primID,primID);
        store8f(valid,&ray.Ng.x,Ng.x);
        store8f(valid,&ray.Ng.y,Ng.y);
        store8f(valid,&ray.Ng.z,Ng.z);
      }
    }

    static __forceinline void intersect(const avxb& valid, Ray8& ray, const Triangle4i* tri, size_t num, const void* geom)
    {
      for (size_t i=0; i<num; i++)
        intersect(valid,ray,tri[i],geom);
    }
    
    static __forceinline avxb occluded(const avxb& valid_i, Ray8& ray, const Triangle4i& tri, const void* geom)
    {
      avxb valid0 = valid_i;

      for (size_t i=0; i<tri.size(); i++)
      {
        STAT3(shadow.trav_prims,1,popcnt(valid_i),8);

        /* load vertices */
        const Vec3f& p0 = *tri.v0[i];
        const Vec3f& p1 = *(Vec3f*)((int*)&p0 + tri.v1[i]);
        const Vec3f& p2 = *(Vec3f*)((int*)&p0 + tri.v2[i]);

        /* calculate vertices relative to ray origin */
        avxb valid = valid0;
        const avx3f O = ray.org;
        const avx3f D = ray.dir;
        const avx3f v0 = avx3f(p0)-O;
        const avx3f v1 = avx3f(p1)-O;
        const avx3f v2 = avx3f(p2)-O;

        /* calculate triangle edges */
        const avx3f e0 = v2-v0;
        const avx3f e1 = v0-v1;
        const avx3f e2 = v1-v2;
        
        /* calculate geometry normal and denominator */
        const avx3f Ng1 = cross(e1,e0);
        const avx3f Ng = Ng1+Ng1;
        const avxf den = dot(avx3f(Ng),D);
        const avxf absDen = abs(den);
        const avxf sgnDen = signmsk(den);
        
        /* perform edge tests */
        const avxf U = dot(avx3f(cross(v2+v0,e0)),D) ^ sgnDen;
        valid &= U >= 0.0f;
        if (likely(none(valid))) continue;
        const avxf V = dot(avx3f(cross(v0+v1,e1)),D) ^ sgnDen;
        valid &= V >= 0.0f;
        if (likely(none(valid))) continue;
        const avxf W = dot(avx3f(cross(v1+v2,e2)),D) ^ sgnDen;
        valid &= W >= 0.0f;
        if (likely(none(valid))) continue;
        
        /* perform depth test */
        const avxf T = dot(v0,avx3f(Ng)) ^ sgnDen;
        valid &= (T >= absDen*ray.tnear) & (absDen*ray.tfar >= T);

        /* perform backface culling */
#if defined(__BACKFACE_CULLING__)
        valid &= den > avxf(zero);
        if (unlikely(none(valid))) continue;
#else
        valid &= den != avxf(zero);
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
        if (unlikely(geometry->hasOcclusionFilter8()))
        {
          /* calculate hit information */
          const avxf u = U / absDen;
          const avxf v = V / absDen;
          const avxf t = T / absDen;
          const int primID = tri.primID[i];
          valid = runOcclusionFilter8(valid,geometry,ray,u,v,t,Ng,geomID,primID);
        }
#endif

        /* update occlusion */
        valid0 &= !valid;
        if (none(valid0)) break;
      }
      return !valid0;
    }

    static __forceinline avxb occluded(const avxb& valid, Ray8& ray, const Triangle4i* tri, size_t num, const void* geom)
    {
      avxb valid0 = valid;
      for (size_t i=0; i<num; i++) {
        valid0 &= !occluded(valid0,ray,tri[i],geom);
        if (none(valid0)) break;
      }
      return !valid0;
    }
  };
}

#endif


