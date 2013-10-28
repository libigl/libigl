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

#ifndef __EMBREE_INSTANCE_H__
#define __EMBREE_INSTANCE_H__

#include "handle.h"

#include "../shapes/shape.h"
#include "../materials/material.h"
#include "../lights/light.h"

namespace embree
{
  /*! Primitive Handle */
  class PrimitiveHandle : public _RTHandle {
    ALIGNED_CLASS;
  public:
    
    /*! Constructs new shape primitive. */
    PrimitiveHandle (const Ref<InstanceHandle<Shape> >& shape, 
                     const Ref<InstanceHandle<Material> >& material, 
                     const AffineSpace3f& transform)
      : bstatic(false), shape(shape), light(null), material(material), transform(transform),
        illumMask(-1), shadowMask(-1) { }

    /*! Constructs new light primitive. */
    PrimitiveHandle (const Ref<InstanceHandle<Light> >& light, 
                     const Ref<InstanceHandle<Material> >& material, 
                     const AffineSpace3f& transform)
      : bstatic(false), shape(null), light(light), material(material), transform(transform),
        illumMask(-1), shadowMask(-1) { }

    /*! Constructs new primitive. */
    PrimitiveHandle (const AffineSpace3f& transform, const Ref<PrimitiveHandle>& other)
      : bstatic(other->bstatic), shape(other->shape), light(other->light), material(other->material), 
        transform(transform * other->transform),
        illumMask(other->illumMask), shadowMask(other->shadowMask) { }

    /*! Setting parameters. */
    void set(const std::string& property, const Variant& data) 
    { 
      if      (property == "illumMask" ) { illumMask  = data.getInt(); }
      else if (property == "shadowMask") { shadowMask = data.getInt(); }
      else if (property == "static"    ) { bstatic    = data.getBool(); }
    }

    void create () {}
    
    Ref<Shape> getShapeInstance() { if (!shape) return NULL; return shape->getInstance(); }
    Ref<Light> getLightInstance() { if (!light) return NULL; return light->getInstance(); }
    Ref<Material> getMaterialInstance() { if (!material) return NULL; return material->getInstance(); }

    /*! input */
  public: 
    bool bstatic;
    Ref<InstanceHandle<Shape> > shape;          //!< Shape in case of a shape primitive
    Ref<InstanceHandle<Light> > light;       //!< Light in case of a light primitive
    Ref<InstanceHandle<Material> > material; //!< Material of shape primitive
    AffineSpace3f transform;             //!< Transformation of primitive
    light_mask_t illumMask;
    light_mask_t shadowMask;
  };
}

#endif
