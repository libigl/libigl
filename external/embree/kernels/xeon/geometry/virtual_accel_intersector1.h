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

#include "virtual_accel.h"
#include "common/ray.h"

namespace embree
{
  namespace isa
  {
    struct VirtualAccelIntersector1
    {
      typedef AccelSetItem Primitive;
      
      struct Precalculations {
        __forceinline Precalculations (const Ray& ray) {}
      };
      
      static __forceinline void intersect(const Precalculations& pre, Ray& ray, const Primitive& prim, Scene* scene) 
      {
        AVX_ZERO_UPPER();
        prim.accel->intersect((RTCRay&)ray,prim.item);
      }
      
      static __forceinline bool occluded(const Precalculations& pre, Ray& ray, const Primitive& prim, Scene* scene) 
      {
        AVX_ZERO_UPPER();
        prim.accel->occluded((RTCRay&)ray,prim.item);
        return ray.geomID == 0;
      }
    };
  }
}
