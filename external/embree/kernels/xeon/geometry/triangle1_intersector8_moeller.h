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

#include "triangle1.h"
#include "common/ray8.h"
#include "geometry/filter.h"

namespace embree
{
  namespace isa
  {
    /*! Intersector for individual precomputed triangles with 8
     *  rays. This intersector implements a modified version of the
     *  Moeller Trumbore intersector from the paper "Fast, Minimum
     *  Storage Ray-Triangle Intersection". In contrast to the paper we
     *  precalculate some factors and factor the calculations
     *  differently to allow precalculating the cross product e1 x
     *  e2. */
    template<bool list>
      struct Triangle1Intersector8MoellerTrumbore
      {
        typedef Triangle1 Primitive;
        
        struct Precalculations {
          __forceinline Precalculations (const avxb& valid, const Ray8& ray) {}
        };
        
        static __forceinline void intersect(const avxb& valid_i, Precalculations& pre, Ray8& ray, const Primitive& tri, Scene* scene)
        {
          STAT3(normal.trav_prims,1,popcnt(valid_i),8);
          
          avxb valid = valid_i;
          const avx3f org = ray.org;
          const avx3f dir = ray.dir;
          
          /* load vertices and calculate edges */
          const avxf v0 = broadcast4f(&tri.v0);
          const avxf v1 = broadcast4f(&tri.v1);
          const avxf v2 = broadcast4f(&tri.v2);
          const avxf e1 = v0-v1;
          const avxf e2 = v2-v0;
          
          /* calculate denominator */
          const avx3f _v0 = avx3f(shuffle<0>(v0),shuffle<1>(v0),shuffle<2>(v0));
          const avx3f C =  _v0 - org;
          const avx3f Ng = avx3f(tri.Ng);
          const avxf den = dot(dir,Ng);
          const avxf sgnDen = signmsk(den);
          const avxf absDen = abs(den);
#if defined(RTCORE_BACKFACE_CULLING)
          valid &= den > avxf(zero);
#else
          valid &= den != avxf(zero);
#endif
          
          /* perform edge tests */
          const avx3f R = cross(dir,C);
          const avx3f _e1(shuffle<0>(e1),shuffle<1>(e1),shuffle<2>(e1));
          const avxf V = dot(R,_e1)^sgnDen;
          const avx3f _e2(shuffle<0>(e2),shuffle<1>(e2),shuffle<2>(e2));
          const avxf U = dot(R,_e2)^sgnDen;
          valid &= V >= avxf(zero) & U >= avxf(zero) & U+V <= absDen;
          if (unlikely(none(valid))) return;
          
          /* perform depth test */
          const avxf T = dot(C,Ng) ^ sgnDen;
          valid &= (T >= absDen*ray.tnear) & (absDen*ray.tfar >= T);
          if (unlikely(none(valid))) return;
          
          /* ray masking test */
#if defined(RTCORE_RAY_MASK)
          valid &= (tri.mask() & ray.mask) != 0;
          if (unlikely(none(valid))) return;
#endif
          
          /* calculate hit information */
          const avxf rcpAbsDen = rcp(absDen);
          const avxf u = U*rcpAbsDen;
          const avxf v = V*rcpAbsDen;
          const avxf t = T*rcpAbsDen;
          const int geomID = tri.geomID<list>();
          const int primID = tri.primID<list>();
          
          /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER)
          Geometry* geometry = scene->get(geomID);
          if (unlikely(geometry->hasIntersectionFilter8())) {
            runIntersectionFilter8(valid,geometry,ray,u,v,t,Ng,geomID,primID);
            return;
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
        
        static __forceinline avxb occluded(const avxb& valid_i, Precalculations& pre, Ray8& ray, const Primitive& tri, Scene* scene)
        {
          STAT3(shadow.trav_prims,1,popcnt(valid_i),8);
          
          avxb valid = valid_i;
          const avx3f org = ray.org;
          const avx3f dir = ray.dir;
          
          /* load vertices and calculate edges */
          const avxf v0 = broadcast4f(&tri.v0);
          const avxf v1 = broadcast4f(&tri.v1);
          const avxf v2 = broadcast4f(&tri.v2);
          const avxf e1 = v0-v1;
          const avxf e2 = v2-v0;
          
          /* calculate denominator */
          const avx3f _v0 = avx3f(shuffle<0>(v0),shuffle<1>(v0),shuffle<2>(v0));
          const avx3f C =  _v0 - org;
          const avxf Ng = broadcast4f(&tri.Ng);
          const avx3f _Ng = avx3f(shuffle<0>(Ng),shuffle<1>(Ng),shuffle<2>(Ng));
          const avxf den = dot(dir,_Ng);
          const avxf sgnDen = signmsk(den);
          const avxf absDen = abs(den);
#if defined(RTCORE_BACKFACE_CULLING)
          valid &= den > avxf(zero);
#else
          valid &= den != avxf(zero);
#endif
          
          /* perform edge tests */
          const avx3f R = cross(dir,C);
          const avx3f _e1(shuffle<0>(e1),shuffle<1>(e1),shuffle<2>(e1));
          const avxf V = dot(R,_e1)^sgnDen;
          const avx3f _e2(shuffle<0>(e2),shuffle<1>(e2),shuffle<2>(e2));
          const avxf U = dot(R,_e2)^sgnDen;
          valid &= V >= avxf(zero) & U >= avxf(zero) & U+V <= absDen;
          if (unlikely(none(valid))) return valid;
          
          /* perform depth test */
          const avxf T = dot(C,_Ng) ^ sgnDen;
          valid &= (T >= absDen*ray.tnear) & (absDen*ray.tfar >= T);
          if (unlikely(none(valid))) return valid;
          
          /* ray masking test */
#if defined(RTCORE_RAY_MASK)
          valid &= (tri.mask() & ray.mask) != 0;
          if (unlikely(none(valid))) return valid;
#endif
          
          /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER)
          const int geomID = tri.geomID<list>();
          Geometry* geometry = scene->get(geomID);
          if (unlikely(geometry->hasOcclusionFilter8()))
          {
            /* calculate hit information */
            const avxf rcpAbsDen = rcp(absDen);
            const avxf u = U*rcpAbsDen;
            const avxf v = V*rcpAbsDen;
            const avxf t = T*rcpAbsDen;
            const int primID = tri.primID<list>();
            valid = runOcclusionFilter8(valid,geometry,ray,u,v,t,Ng,geomID,primID);
          }
#endif
          return valid;
        }
      };
  }
}
