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

#ifndef __EMBREE_LIGHT_H__
#define __EMBREE_LIGHT_H__

#include "api/parms.h"
#include "shapes/shape.h"
#include "shapes/differentialgeometry.h"

namespace embree
{
  /*! Interface to different lights. A light can be evaluated,
   *  sampled, and the sampling PDF be evaluated. */
  class Light : public RefCount {
    ALIGNED_CLASS
  public:

    /*! Lights need a virtual destructor. */
    virtual ~Light() {}

    /*! Instantiates a new light by transforming this light to a different location. */
    virtual Ref<Light> transform (const AffineSpace3f& xfm) const = 0;

    /*! Returns the shape of the light. */
    virtual Ref<Shape> shape() { return null; }

    /*! Evaluates the radiance that would arrive at a given location
     *  from a given direction assuming no blocking geometry. \returns
     *  the emitted radiance. */
    virtual Col3f eval (const DifferentialGeometry& dg, /*!< The shade point that is illuminated.    */
                        const Vec3f& wi)                /*!< The direction the light is coming from. */ const
    {
      return zero;
    }

    /*! Samples the light for a point to shade. \returns the radiance
     *  arriving from the sampled direction. */
    virtual Col3f sample (const DifferentialGeometry& dg, /*!< The shade point that is illuminated. */
                          Sample3f& wi,                   /*!< Returns the sampled direction including PDF.*/
                          float& tMax,                    /*!< Returns the distance of the light. */
                          const Vec2f& sample)            /*!< Sample locations are provided by the caller. */ const = 0;

    /*! Evaluates the probability distribution function used by the
     *  sampling function of the light for a shade location and
     *  direction. \returns the probability density */
    virtual float pdf (const DifferentialGeometry& dg,   /*!< The shade location to compute the PDF for. */
                       const Vec3f& wi)                   /*!< The direction to compute the PDF for. */ const = 0;

    /*! Indicates that the sampling of the light is expensive and the
     *  integrator should presample the light. */
    virtual bool precompute() const { return false; }
  };

  /*! Interface to an area light. In addition to a basic light, the
   * area light allows to query the emitted radiance at a location on
   * the light itself. */
  class AreaLight : public Light {
  public:

    /*! Returns the the emitted radiance for a location on the light
     *  and direction. */
    virtual Col3f Le(const DifferentialGeometry& dg,     /*!< The location on the light. */
                     const Vec3f& wo)                    /*!< The direction the light is emitted into. */ const = 0;
  };

  /*! Interface to an environment light. In addition to a basic light,
   * an environment light allows to query the emitted radiance for a
   * direction. */
  class EnvironmentLight : public Light {
  public:

    /*! Returns the emitted radiance of the environment light. */
    virtual Col3f Le(const Vec3f& wo                     /*!< The direction the light comes from. */) const = 0;
  };
}

#endif
