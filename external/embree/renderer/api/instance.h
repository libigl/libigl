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

#ifndef __EMBREE_INSTANCE_H__
#define __EMBREE_INSTANCE_H__

#include "shapes/shape.h"
#include "materials/material.h"
#include "lights/light.h"

namespace embree
{
  /*! Shape instance. A shape instance is a shape with attached
   *  material and area light. */
  class Instance : public RefCount
  {
  public:

    /*! Shape instance constructor. */
    Instance (size_t id,                       /*!< User ID of the instance.          */
              const Ref<Shape>& shape,         /*!< Shape of the instance.            */
              const Ref<Material>& material,   /*!< Material attached to the shape.   */
              const Ref<AreaLight>& light)     /*!< Area light attached to the shape. */
      : id(id), shape(shape), material(material), light(light) {}

    /*! Post intersection for the instance. The instance is
     *  responsible for setting the material and area light, while the
     *  shape for interpolating shading data. */
    __forceinline void postIntersect(const Ray& ray, DifferentialGeometry& dg) const {
      dg.material = material.ptr;
      dg.light = light.ptr;
      shape->postIntersect(ray,dg);
    }

  public:
    size_t id;              //!< User ID of the instance.
    Ref<Shape> shape;       //!< Shape of the instance.
    Ref<Material> material; //!< Material attached to the shape.
    Ref<AreaLight> light;   //!< Area light attached to the shape.
  };
}

#endif
