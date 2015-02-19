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

#include "bezier1i.h"
#include "bezier_intersector8.h"

namespace embree
{
  namespace isa
  {
    /*! Intersector for a single ray from a ray packet with a bezier curve. */
    template<bool list>
      struct Bezier1iIntersector8
      {
        typedef Bezier1i Primitive;
        typedef BezierIntersector8::Precalculations Precalculations;
        
        static __forceinline void intersect(Precalculations& pre, Ray8& ray, const size_t k, const Primitive& curve, Scene* scene) 
        {
          const BezierCurves* in = (BezierCurves*) scene->get(curve.geomID<list>());
          const Vec3fa a0 = in->vertex(curve.vertexID+0,0);
          const Vec3fa a1 = in->vertex(curve.vertexID+1,0);
          const Vec3fa a2 = in->vertex(curve.vertexID+2,0);
          const Vec3fa a3 = in->vertex(curve.vertexID+3,0);
          BezierIntersector8::intersect(pre,ray,k,a0,a1,a2,a3,curve.geomID<list>(),curve.primID<list>(),scene);
        }
        
        static __forceinline void intersect(const avxb& valid_i, Precalculations& pre, Ray8& ray, const Primitive& curve, Scene* scene)
        {
          int mask = movemask(valid_i);
          while (mask) intersect(pre,ray,__bscf(mask),curve,scene);
        }
        
        static __forceinline bool occluded(Precalculations& pre, Ray8& ray, const size_t k, const Primitive& curve, Scene* scene) 
        {
          const BezierCurves* in = (BezierCurves*) scene->get(curve.geomID<list>());
          const Vec3fa a0 = in->vertex(curve.vertexID+0,0);
          const Vec3fa a1 = in->vertex(curve.vertexID+1,0);
          const Vec3fa a2 = in->vertex(curve.vertexID+2,0);
          const Vec3fa a3 = in->vertex(curve.vertexID+3,0);
          return BezierIntersector8::occluded(pre,ray,k,a0,a1,a2,a3,curve.geomID<list>(),curve.primID<list>(),scene);
        }
        
        static __forceinline avxb occluded(const avxb& valid_i, Precalculations& pre, Ray8& ray, const Primitive& curve, Scene* scene)
        {
          avxb valid_o = false;
          int mask = movemask(valid_i);
          while (mask) {
            size_t k = __bscf(mask);
            if (occluded(pre,ray,k,curve,scene))
              valid_o[k] = -1;
          }
          return valid_o;
        }
      };
    
    template<bool list>
      struct Bezier1iIntersector8MB
      {
        typedef Bezier1iMB Primitive;
        typedef BezierIntersector8::Precalculations Precalculations;
        
        static __forceinline void intersect(Precalculations& pre, Ray8& ray, const size_t k, const Primitive& curve, Scene* scene)
        {
          const BezierCurves* in = (BezierCurves*) scene->get(curve.geomID<list>());
          const Vec3fa a0 = in->vertex(curve.vertexID+0,0);
          const Vec3fa a1 = in->vertex(curve.vertexID+1,0);
          const Vec3fa a2 = in->vertex(curve.vertexID+2,0);
          const Vec3fa a3 = in->vertex(curve.vertexID+3,0);
          const Vec3fa b0 = in->vertex(curve.vertexID+0,1);
          const Vec3fa b1 = in->vertex(curve.vertexID+1,1);
          const Vec3fa b2 = in->vertex(curve.vertexID+2,1);
          const Vec3fa b3 = in->vertex(curve.vertexID+3,1);
          const float t0 = 1.0f-ray.time[k], t1 = ray.time[k];
          const Vec3fa p0 = t0*a0 + t1*b0;
          const Vec3fa p1 = t0*a1 + t1*b1;
          const Vec3fa p2 = t0*a2 + t1*b2;
          const Vec3fa p3 = t0*a3 + t1*b3;
          BezierIntersector8::intersect(pre,ray,k,p0,p1,p2,p3,curve.geomID<list>(),curve.primID<list>(),scene);
        }
        
        static __forceinline bool occluded(Precalculations& pre, Ray8& ray, const size_t k, const Primitive& curve, Scene* scene) 
        {
          const BezierCurves* in = (BezierCurves*) scene->get(curve.geomID<list>());
          const Vec3fa a0 = in->vertex(curve.vertexID+0,0);
          const Vec3fa a1 = in->vertex(curve.vertexID+1,0);
          const Vec3fa a2 = in->vertex(curve.vertexID+2,0);
          const Vec3fa a3 = in->vertex(curve.vertexID+3,0);
          const Vec3fa b0 = in->vertex(curve.vertexID+0,1);
          const Vec3fa b1 = in->vertex(curve.vertexID+1,1);
          const Vec3fa b2 = in->vertex(curve.vertexID+2,1);
          const Vec3fa b3 = in->vertex(curve.vertexID+3,1);
          const float t0 = 1.0f-ray.time[k], t1 = ray.time[k];
          const Vec3fa p0 = t0*a0 + t1*b0;
          const Vec3fa p1 = t0*a1 + t1*b1;
          const Vec3fa p2 = t0*a2 + t1*b2;
          const Vec3fa p3 = t0*a3 + t1*b3;
          return BezierIntersector8::occluded(pre,ray,k,p0,p1,p2,p3,curve.geomID<list>(),curve.primID<list>(),scene);
        }
      };
  }
}
