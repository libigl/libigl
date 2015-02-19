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

#include "triangle4i.h"
#include "common/ray.h"
#include "common/scene_triangle_mesh.h"

namespace embree
{
  namespace isa
  {
    /*! Intersector1 for triangle4i */
    template<bool list>
      struct Triangle4iIntersector1Pluecker
      {
        typedef Triangle4i Primitive;
        
        struct Precalculations {
          __forceinline Precalculations (const Ray& ray) {}
        };
        
        static __forceinline void intersect(const Precalculations& pre, Ray& ray, const Primitive& tri, Scene* scene)
        {
          /* gather vertices */
          STAT3(normal.trav_prims,1,1,1);
          const int* base0 = (const int*) tri.v0[0];
          const int* base1 = (const int*) tri.v0[1];
          const int* base2 = (const int*) tri.v0[2];
          const int* base3 = (const int*) tri.v0[3];
          const ssef a0 = loadu4f(base0          ), a1 = loadu4f(base1          ), a2 = loadu4f(base2          ), a3 = loadu4f(base3          );
          const ssef b0 = loadu4f(base0+tri.v1[0]), b1 = loadu4f(base1+tri.v1[1]), b2 = loadu4f(base2+tri.v1[2]), b3 = loadu4f(base3+tri.v1[3]);
          const ssef c0 = loadu4f(base0+tri.v2[0]), c1 = loadu4f(base1+tri.v2[1]), c2 = loadu4f(base2+tri.v2[2]), c3 = loadu4f(base3+tri.v2[3]);
          sse3f p0; transpose(a0,a1,a2,a3,p0.x,p0.y,p0.z);
          sse3f p1; transpose(b0,b1,b2,b3,p1.x,p1.y,p1.z);
          sse3f p2; transpose(c0,c1,c2,c3,p2.x,p2.y,p2.z);
          
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
          
          /* calculate geometry normal and denominator */
          const sse3f Ng1 = cross(e1,e0);
          const sse3f Ng = Ng1+Ng1;
          const ssef den = dot(Ng,D);
          const ssef absDen = abs(den);
          const ssef sgnDen = signmsk(den);
          
          /* perform edge tests */
          const ssef U = dot(cross(v2+v0,e0),D) ^ sgnDen;
          const ssef V = dot(cross(v0+v1,e1),D) ^ sgnDen;
          const ssef W = dot(cross(v1+v2,e2),D) ^ sgnDen;
          sseb valid = (U >= 0.0f) & (V >= 0.0f) & (W >= 0.0f);
          if (unlikely(none(valid))) return;
          
          /* perform depth test */
          const ssef T = dot(v0,Ng) ^ sgnDen;
          valid &= (T >= absDen*ssef(ray.tnear)) & (absDen*ssef(ray.tfar) >= T);
          if (unlikely(none(valid))) return;
          
          /* perform backface culling */
#if defined(RTCORE_BACKFACE_CULLING)
          valid &= den > ssef(zero);
          if (unlikely(none(valid))) return;
#else
          valid &= den != ssef(zero);
          if (unlikely(none(valid))) return;
#endif
          
          /* calculate hit information */
          const ssef u = U / absDen;
          const ssef v = V / absDen;
          const ssef t = T / absDen;
          size_t i = select_min(valid,t);
          int geomID = tri.geomID<list>(i);
          
          /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER) || defined(RTCORE_RAY_MASK)
          while (true) 
          {
            TriangleMesh* geometry = scene->getTriangleMesh(geomID);
            if ((geometry->mask & ray.mask) == 0) {
              valid[i] = 0;
              if (none(valid)) return;
              i = select_min(valid,t);
              geomID = tri.geomID<list>(i);
              continue;
            }
            if (likely(!geometry->hasIntersectionFilter1())) 
            {
#endif
              /* update hit information */
              ray.u = u[i];
              ray.v = v[i];
              ray.tfar = t[i];
              ray.Ng.x = Ng.x[i];
              ray.Ng.y = Ng.y[i];
              ray.Ng.z = Ng.z[i];
              ray.geomID = geomID;
              ray.primID = tri.primID<list>(i);
              
#if defined(RTCORE_INTERSECTION_FILTER) || defined(RTCORE_RAY_MASK)
              return;
            }
            
            Vec3fa N = Vec3fa(Ng.x[i],Ng.y[i],Ng.z[i]);
            if (runIntersectionFilter1(geometry,ray,u[i],v[i],t[i],N,geomID,tri.primID<list>(i))) return;
            valid[i] = 0;
            if (none(valid)) return;
            i = select_min(valid,t);
            geomID = tri.geomID<list>(i);
          }
#endif
          
        }
        
        static __forceinline bool occluded(const Precalculations& pre, Ray& ray, const Primitive& tri, Scene* scene)
        {
          /* gather vertices */
          STAT3(shadow.trav_prims,1,1,1);
          const int* base0 = (const int*) tri.v0[0];
          const int* base1 = (const int*) tri.v0[1];
          const int* base2 = (const int*) tri.v0[2];
          const int* base3 = (const int*) tri.v0[3];
          const ssef a0 = loadu4f(base0          ), a1 = loadu4f(base1          ), a2 = loadu4f(base2          ), a3 = loadu4f(base3          );
          const ssef b0 = loadu4f(base0+tri.v1[0]), b1 = loadu4f(base1+tri.v1[1]), b2 = loadu4f(base2+tri.v1[2]), b3 = loadu4f(base3+tri.v1[3]);
          const ssef c0 = loadu4f(base0+tri.v2[0]), c1 = loadu4f(base1+tri.v2[1]), c2 = loadu4f(base2+tri.v2[2]), c3 = loadu4f(base3+tri.v2[3]);
          sse3f p0; transpose(a0,a1,a2,a3,p0.x,p0.y,p0.z);
          sse3f p1; transpose(b0,b1,b2,b3,p1.x,p1.y,p1.z);
          sse3f p2; transpose(c0,c1,c2,c3,p2.x,p2.y,p2.z);
          
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
          
          /* calculate geometry normal and denominator */
          const sse3f Ng1 = cross(e1,e0);
          const sse3f Ng = Ng1+Ng1;
          const ssef den = dot(Ng,D);
          const ssef absDen = abs(den);
          const ssef sgnDen = signmsk(den);
          
          /* perform edge tests */
          const ssef U = dot(cross(v2+v0,e0),D) ^ sgnDen;
          const ssef V = dot(cross(v0+v1,e1),D) ^ sgnDen;
          const ssef W = dot(cross(v1+v2,e2),D) ^ sgnDen;
          sseb valid = (U >= 0.0f) & (V >= 0.0f) & (W >= 0.0f);
          if (unlikely(none(valid))) return false;
          
          /* perform depth test */
          const ssef T = dot(v0,Ng) ^ sgnDen;
          valid &= (T >= absDen*ssef(ray.tnear)) & (absDen*ssef(ray.tfar) >= T);
          if (unlikely(none(valid))) return false;
          
          /* perform backface culling */
#if defined(RTCORE_BACKFACE_CULLING)
          valid &= den > ssef(zero);
          if (unlikely(none(valid))) return false;
#else
          valid &= den != ssef(zero);
          if (unlikely(none(valid))) return false;
#endif
          
          /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER) || defined(RTCORE_RAY_MASK)
          
          for (size_t m=movemask(valid), i=__bsf(m); m!=0; m=__btc(m,i), i=__bsf(m))
          {  
            const int geomID = tri.geomID<list>(i);
            TriangleMesh* geometry = scene->getTriangleMesh(geomID);
            
#if defined(RTCORE_RAY_MASK)
            /* goto next hit if mask test fails */
            if ((geometry->mask & ray.mask) == 0) 
              continue;
#endif
            
#if defined(RTCORE_INTERSECTION_FILTER)
            
            /* if we have no filter then the test passes */
            if (likely(!geometry->hasOcclusionFilter1()))
              return true;
            
            /* calculate hit information */
            const ssef rcpAbsDen = rcp(absDen);
            const ssef u = U * rcpAbsDen;
            const ssef v = V * rcpAbsDen;
            const ssef t = T * rcpAbsDen;
            const Vec3fa N = Vec3fa(Ng.x[i],Ng.y[i],Ng.z[i]);
            if (runOcclusionFilter1(geometry,ray,u[i],v[i],t[i],N,geomID,tri.primID<list>(i))) 
#endif
              return true;
          }
          return false;
#else
          return true;
#endif
        }
      };
  }
}
