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

#include "bezier1v.h"
#include "bezier_intersector4.h"

namespace embree
{
  namespace isa
  {
    /*! Intersector for a single ray from a ray packet with a bezier curve. */
    template<bool list>
      struct Bezier1vIntersector4
      {
        typedef Bezier1v Primitive;
        typedef BezierIntersector4::Precalculations Precalculations;
        
        static __forceinline void intersect(Precalculations& pre, Ray4& ray, const size_t k, const Primitive& curve, Scene* scene) {
          BezierIntersector4::intersect(pre,ray,k,curve.p0,curve.p1,curve.p2,curve.p3,curve.geomID<list>(),curve.primID<list>(),scene);
        }
        
        static __forceinline void intersect(const sseb& valid_i, Precalculations& pre, Ray4& ray, const Primitive& curve, Scene* scene)
        {
          int mask = movemask(valid_i);
          while (mask) intersect(pre,ray,__bscf(mask),curve,scene);
        }
        
        static __forceinline bool occluded(Precalculations& pre, Ray4& ray, const size_t k, const Primitive& curve, Scene* scene) {
          return BezierIntersector4::occluded(pre,ray,k,curve.p0,curve.p1,curve.p2,curve.p3,curve.geomID<list>(),curve.primID<list>(),scene);
        }
        
        static __forceinline sseb occluded(const sseb& valid_i, Precalculations& pre, Ray4& ray, const Primitive& curve, Scene* scene)
        {
          sseb valid_o = false;
          int mask = movemask(valid_i);
          while (mask) {
            size_t k = __bscf(mask);
            if (occluded(pre,ray,k,curve,scene))
              valid_o[k] = -1;
          }
          return valid_o;
        }
      };
  }
}

