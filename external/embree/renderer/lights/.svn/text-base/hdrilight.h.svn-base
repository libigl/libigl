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

#ifndef __EMBREE_HDRI_LIGHT_H__
#define __EMBREE_HDRI_LIGHT_H__

#include "lights/light.h"
#include "image/image.h"
#include "samplers/distribution2d.h"

namespace embree
{
  /*! Implements a texture mapped environment light. */
  class HDRILight : public EnvironmentLight
  {
    /*! construction from members */
    HDRILight(const AffineSpace3f& local2world, unsigned width, unsigned height, 
              const Col3f& L, const Ref<Image>& pixels, const Ref<Distribution2D>& distribution)
      : local2world(local2world), world2local(rcp(local2world)), width(width), height(height), 
        L(L), pixels(pixels), distribution(distribution) {}

  public:

    /*! Construction from parameter container. */
    HDRILight (const Parms& parms);

    Ref<Light> transform(const AffineSpace3f& xfm) const {
      return new HDRILight(xfm * local2world,width,height,L,pixels,distribution);
    }

    Col3f Le    (const Vec3f& wo) const;
    Col3f eval  (const DifferentialGeometry& dg, const Vec3f& wi) const;
    Col3f sample(const DifferentialGeometry& dg, Sample3f& wi, float& tMax, const Vec2f& s) const;
    float pdf   (const DifferentialGeometry& dg, const Vec3f& wi) const;
    bool precompute() const { return true; }

  private:
    AffineSpace3f local2world;            //!< Transformation from light space into world space
    AffineSpace3f world2local;            //!< Transformation from world space into light space
    unsigned width, height;             //!< Width and height of the used image.
    Col3f L;                            //!< Scaling factor for the image.
    Ref<Image> pixels;                  //!< The image mapped to the environment.
    Ref<Distribution2D> distribution;   //!< The 2D distribution used to importance sample the image.
  };
}

#endif
