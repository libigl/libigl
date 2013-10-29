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

#ifndef __EMBREE_ACCEL_VIRTUAL_OBJECT_INTERSECTOR1_H__
#define __EMBREE_ACCEL_VIRTUAL_OBJECT_INTERSECTOR1_H__

#include "virtual_object.h"
#include "../common/ray.h"

namespace embree
{
  struct VirtualObjectIntersector1
  {
    typedef VirtualObject Triangle;

    /*! Name of intersector */
    static const char* name() { return "default"; }
    
    static __forceinline void intersect(Ray& ray, const VirtualObject& id, const Vec3fa* vertices) 
    {
      const VirtualScene::Object& obj = ((const VirtualScene::Object*) vertices)[id.id];

      /*! skip this object if masked out by the ray */
      if ((obj.mask & ray.mask) == 0) return;

      /*! fast path for identity transformation */
      if (likely(obj.hasTransform == false)) {
        obj.intersector1->intersect(ray);
        if (obj.id0 != 0x7FFFFFFF && ray.id0 == 0x7FFFFFFF) 
          ray.id0 = obj.id0; 
        return;
      }

      /*! slow path with full transformation */
      const Vector3f org = ray.org, dir = ray.dir;
      ray.org = xfmPoint (obj.world2local,org);
      ray.dir = xfmVector(obj.world2local,dir);
      obj.intersector1->intersect(ray);
      ray.org = org;
      ray.dir = dir;

      if (obj.id0 != 0x7FFFFFFF && ray.id0 == 0x7FFFFFFF)
        ray.id0 = obj.id0; 
    }

    static __forceinline void intersect(Ray& ray, const VirtualObject* id, size_t items, const Vec3fa* vertices) 
    {
      for (size_t i=0; i<items; i++)
        intersect(ray,id[i],vertices);
    }

    static __forceinline bool occluded(Ray& ray, const VirtualObject& id, const Vec3fa* vertices) 
    {
      const VirtualScene::Object& obj = ((const VirtualScene::Object*) vertices)[id.id];
      
      /*! skip this object if masked out by the ray */
      if ((obj.mask & ray.mask) == 0) return false;

      /*! fast path for identity transformation */
      if (likely(obj.hasTransform == false))
        return obj.intersector1->occluded(ray);

      /*! slow path with full transformation */
      const Vector3f org = ray.org, dir = ray.dir;
      ray.org = xfmPoint (obj.world2local,org);
      ray.dir = xfmVector(obj.world2local,dir);
      bool ret = obj.intersector1->occluded(ray);
      ray.org = org;
      ray.dir = dir;
      return ret;
    }

    static __forceinline bool occluded(Ray& ray, const VirtualObject* id, size_t items, const Vec3fa* vertices) 
    {
      for (size_t i=0; i<items; i++)
        if (occluded(ray,id[i],vertices))
          return true;
      return false;
    }
  };
}

#endif


