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

#include "common/accel.h"
#include "common/ray.h"

namespace embree
{
  struct VirtualAccelIntersector1
  {
    typedef AccelSetItem Primitive;

    static __forceinline void intersect(Ray& ray, const Primitive& prim, const void* geom) 
    {
      prim.accel->intersect((RTCRay&)ray,prim.item);
    }

    static __forceinline void intersect(Ray& ray, const Primitive* prim, size_t num, const void* geom) 
    {
      for (size_t i=0; i<num; i++) 
        intersect(ray,prim[i],geom);
    }

    static __forceinline bool occluded(Ray& ray, const Primitive& prim, const void* geom) 
    {
      prim.accel->occluded((RTCRay&)ray,prim.item);
      return ray.geomID == 0;
    }

    static __forceinline bool occluded(Ray& ray, const Primitive* prim, size_t num, const void* geom) 
    {
      for (size_t i=0; i<num; i++) 
        if (occluded(ray,prim[i],geom))
          return true;

      return false;
    }
  };
}
