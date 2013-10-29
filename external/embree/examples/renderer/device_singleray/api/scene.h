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

#ifndef __EMBREE_BACKEND_SCENE_H__
#define __EMBREE_BACKEND_SCENE_H__

#include "handle.h"
#include "instance.h"

#include "../lights/light.h"
#include "../shapes/differentialgeometry.h"

/*! include interface to ray tracing core */
#include "embree/include/embree.h"
#include "embree/include/intersector1.h"

namespace embree
{
  /*! Scene holding all geometry and lights. */
  class BackendScene : public RefCount
  {
  public:

    class Handle : public InstanceHandle<BackendScene> {
      ALIGNED_CLASS;
    public:
      
      Handle () : accelTy("default"), builderTy("default"), traverserTy("default") {}
      
      void set(const std::string& property, const Variant& data)
      {
        if      (property == "accel") accelTy = data.getString();
        else if (property == "builder") builderTy = data.getString();
        else if (property == "traverser") traverserTy = data.getString();
      }
      
      virtual void setPrimitive(size_t slot, Ref<PrimitiveHandle> prim) = 0;
      
    public:
      std::string accelTy;
      std::string builderTy;
      std::string traverserTy;
    };

  public:

    BackendScene (RTCGeometry* accel, RTCIntersector1* intersector1)
      : accel(accel), intersector(intersector1) {}

    ~BackendScene () {
      if (accel) rtcDeleteGeometry(accel);
      if (intersector) rtcDeleteIntersector1(intersector);
    }

    /*! Adds a light to the scene. */
    void add(const Ref<Light>& light) {
      allLights.push_back(light);
      if (Ref<EnvironmentLight> envlight = dynamic_cast<EnvironmentLight*>(light.ptr)) envLights.push_back(envlight);
    }

    /*! Helper to call the post intersector of the shape instance,
     *  which will call the post intersector of the shape. */
    virtual void postIntersect(const Ray& ray, DifferentialGeometry& dg) const = 0;

  public:
    std::vector<Ref<Light> > allLights;              //!< All lights of the scene
    std::vector<Ref<EnvironmentLight> > envLights;   //!< Environment lights of the scene
    RTCGeometry* accel;                              //!< Acceleration structure over geometry
    RTCIntersector1* intersector;                    //!< Acceleration structure intersector
  };
}

#endif
