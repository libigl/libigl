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

#ifndef __EMBREE_ACCEL_VIRTUAL_OBJECT_INTERSECTOR4_H__
#define __EMBREE_ACCEL_VIRTUAL_OBJECT_INTERSECTOR4_H__

#include "virtual_object.h"
#include "../common/ray4.h"

namespace embree
{
  struct VirtualObjectIntersector4
  {
    typedef VirtualObject Triangle;
  
    /*! Name of intersector */
    static const char* name() { return "default"; }
      
    static __forceinline void intersect(sseb& valid, Ray4& ray, const VirtualObject& id, const Vec3fa* vertices) 
    {
      const VirtualScene::Object& obj = ((const VirtualScene::Object*) vertices)[id.id];

      /*! skip this object if masked out by the ray */
      valid &= (obj.mask & ray.mask) != ssei(zero);
      if (none(valid)) return;

      /*! fast path for identity transformation */
      if (likely(obj.hasTransform == false)) {
        obj.intersector4->intersect(valid,ray);
        if (obj.id0 != 0x7FFFFFFF) 
          ray.id0 = select(ray.id0 == 0x7FFFFFFF, obj.id0, ray.id0);
        return;
      }

      /*! slow path with full transformation */
      const sse3f org = ray.org, dir = ray.dir;
      ray.org = xfmPoint (AffineSpaceT<LinearSpace3<sse3f > >(obj.world2local),org);
      ray.dir = xfmVector(AffineSpaceT<LinearSpace3<sse3f > >(obj.world2local),dir);
      obj.intersector4->intersect(valid,ray);
      ray.org = org;
      ray.dir = dir;

      if (obj.id0 != 0x7FFFFFFF)
        ray.id0 = select(ray.id0 == 0x7FFFFFFF, obj.id0, ray.id0);
    }

    static __forceinline void intersect(sseb& valid, Ray4& ray, const VirtualObject* id, size_t items, const Vec3fa* vertices) 
    {
      for (size_t i=0; i<items; i++)
        intersect(valid,ray,id[i],vertices);
    }

    static __forceinline sseb occluded(sseb& valid, Ray4& ray, const VirtualObject& id, const Vec3fa* vertices) 
    {
      const VirtualScene::Object& obj = ((const VirtualScene::Object*) vertices)[id.id];
      
      /*! skip this object if masked out by the ray */
      valid &= (obj.mask & ray.mask) != ssei(zero);
      if (none(valid)) return false;

      /*! fast path for identity transformation */
      if (likely(obj.hasTransform == false))
        return obj.intersector4->occluded(valid,ray);

      /*! slow path with full transformation */
      const sse3f org = ray.org, dir = ray.dir;
      ray.org = xfmPoint (AffineSpaceT<LinearSpace3<sse3f > >(obj.world2local),org);
      ray.dir = xfmVector(AffineSpaceT<LinearSpace3<sse3f > >(obj.world2local),dir);
      sseb ret = obj.intersector4->occluded(valid,ray);
      ray.org = org;
      ray.dir = dir;
      return ret;
    }

    static __forceinline sseb occluded(sseb& valid, Ray4& ray, const VirtualObject* id, size_t items, const Vec3fa* vertices) 
    {
      sseb not_occluded = valid;
      for (size_t i=0; i<items; i++) {
        not_occluded &= !occluded(not_occluded,ray,id[i],vertices);
        if (none(not_occluded)) return valid;
      }
      return valid & !not_occluded;
    }
  };
}

#endif
