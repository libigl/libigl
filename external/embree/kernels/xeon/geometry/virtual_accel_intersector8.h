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

#ifndef __EMBREE_VIRTUAL_ACCEL_INTERSECTOR8_H__
#define __EMBREE_VIRTUAL_ACCEL_INTERSECTOR8_H__

#include "common/accel.h"
#include "common/ray8.h"

namespace embree
{
  struct VirtualAccelIntersector8
  {
    typedef AccelSetItem Primitive;

    static __forceinline void intersect(const avxb& valid_i, Ray8& ray, const Primitive& prim, const void* geom) 
    {
      AVX_ZERO_UPPER();
      prim.accel->intersect8(&valid_i,(RTCRay8&)ray,prim.item);
    }

    static __forceinline void intersect(const avxb& valid, Ray8& ray, const Primitive* tri, size_t num, const void* geom)
    {
      for (size_t i=0; i<num; i++)
        intersect(valid,ray,tri[i],geom);
    }

    static __forceinline avxb occluded(const avxb& valid_i, const Ray8& ray, const Primitive& prim, const void* geom) 
    {
      AVX_ZERO_UPPER();
      prim.accel->occluded8(&valid_i,(RTCRay8&)ray,prim.item);
      return ray.geomID == 0;
    }

    static __forceinline avxb occluded(const avxb& valid, const Ray8& ray, const Primitive* tri, size_t num, const void* geom)
    {
      avxb terminated = !valid;
      for (size_t i=0; i<num; i++) {
        terminated |= occluded(!terminated,ray,tri[i],geom);
        if (all(terminated)) return terminated;
      }
      return terminated;
    }
  };
}

#endif


