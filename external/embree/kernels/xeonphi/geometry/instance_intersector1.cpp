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

#include "instance_intersector1.h"

namespace embree
{
  namespace isa
  {
    void InstanceBoundsFunction(const Instance* instance, size_t item, BBox3fa& bounds_o)
    {
      Vec3fa lower = instance->object->bounds.lower;
      Vec3fa upper = instance->object->bounds.upper;
      AffineSpace3fa local2world = instance->local2world;
      Vec3fa p000 = xfmPoint(local2world,Vec3fa(lower.x,lower.y,lower.z));
      Vec3fa p001 = xfmPoint(local2world,Vec3fa(lower.x,lower.y,upper.z));
      Vec3fa p010 = xfmPoint(local2world,Vec3fa(lower.x,upper.y,lower.z));
      Vec3fa p011 = xfmPoint(local2world,Vec3fa(lower.x,upper.y,upper.z));
      Vec3fa p100 = xfmPoint(local2world,Vec3fa(upper.x,lower.y,lower.z));
      Vec3fa p101 = xfmPoint(local2world,Vec3fa(upper.x,lower.y,upper.z));
      Vec3fa p110 = xfmPoint(local2world,Vec3fa(upper.x,upper.y,lower.z));
      Vec3fa p111 = xfmPoint(local2world,Vec3fa(upper.x,upper.y,upper.z));
      bounds_o.lower = min(min(min(p000,p001),min(p010,p011)),min(min(p100,p101),min(p110,p111)));
      bounds_o.upper = max(max(max(p000,p001),max(p010,p011)),max(max(p100,p101),max(p110,p111)));
    }

    RTCBoundsFunc InstanceBoundsFunc = (RTCBoundsFunc) InstanceBoundsFunction;

    void FastInstanceIntersector1::intersect(const Instance* instance, Ray& ray, size_t item)
    {
      const Vec3fa ray_org = ray.org;
      const Vec3fa ray_dir = ray.dir;
      const int ray_geomID = ray.geomID;
      const int ray_instID = ray.instID;
      ray.org = xfmPoint (instance->world2local,ray_org);
      ray.dir = xfmVector(instance->world2local,ray_dir);
      ray.geomID = -1;
      ray.instID = instance->id;
      instance->object->intersect((RTCRay&)ray);
      ray.org = ray_org;
      ray.dir = ray_dir;
      if (ray.geomID == -1) {
        ray.geomID = ray_geomID;
        ray.instID = ray_instID;
      }
    }
    
    void FastInstanceIntersector1::occluded (const Instance* instance, Ray& ray, size_t item)
    {
      const Vec3fa ray_org = ray.org;
      const Vec3fa ray_dir = ray.dir;
      ray.org = xfmPoint (instance->world2local,ray_org);
      ray.dir = xfmVector(instance->world2local,ray_dir);
      instance->object->occluded((RTCRay&)ray);
      ray.org = ray_org;
      ray.dir = ray_dir;
    }
    
    DEFINE_SET_INTERSECTOR1(InstanceIntersector1,FastInstanceIntersector1);
  }
}
