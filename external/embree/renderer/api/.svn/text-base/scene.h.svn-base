// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#include "instance.h"
#include "lights/light.h"
#include "shapes/differentialgeometry.h"

/*! include interface to ray tracing core */
#include "rtcore/common/accel.h"
#include "rtcore/common/intersector.h"

namespace embree
{
  /*! Scene holding all geometry and lights. */
  class BackendScene : public RefCount
  {
  public:

    /*! Construction of empty scene. */
    BackendScene () {}

    /*! Adds a light to the scene. */
    void add(const Ref<Light>& light) {
      allLights.push_back(light);
      if (Ref<EnvironmentLight> envlight = dynamic_cast<EnvironmentLight*>(light.ptr)) envLights.push_back(envlight);
    }

    /*! Adds a shape instance to the scene. */
    size_t add(const Ref<Instance>& instance) {
      geometry.push_back(instance);
      return geometry.size()-1;
    }

    /*! Helper to call the post intersector of the shape instance,
     *  which will call the post intersector of the shape. */
    __forceinline void postIntersect(const Ray& ray, DifferentialGeometry& dg) const {
      if (dg) geometry[dg.id0]->postIntersect(ray,dg);
    }

    void setAccel(Ref<Accel> accel) {
      this->accel = accel;
      this->intersector  = accel->queryInterface<Intersector> ();
    }

  public:
    std::vector<Ref<Light> > allLights;              //!< All lights of the scene
    std::vector<Ref<EnvironmentLight> > envLights;   //!< Environemnt lights of the scene
    std::vector<Ref<Instance> > geometry;            //!< Geometries of the scene
    Ref<Accel> accel;                                //!< Acceleration structure over geometry
    Ref<Intersector> intersector;
  };
}

#endif
