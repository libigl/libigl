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

#ifndef __EMBREE_AMBIENT_LIGHT_H__
#define __EMBREE_AMBIENT_LIGHT_H__

#include "../lights/light.h"

namespace embree
{
  /*! Implements an ambient light. An ambient light behaves like a
   *  uniform environment map. */
  class AmbientLight : public EnvironmentLight
  {
  protected:
    /*! Construction from radiance. */
    AmbientLight (const Color& L,
                  light_mask_t illumMask=-1,
                  light_mask_t shadowMask=-1)
      : EnvironmentLight(illumMask,shadowMask), L(L) { }

  public:

    /*! Construction from parameter container. */
    AmbientLight (const Parms& parms) {
      L = parms.getColor("L");
    }

    Ref<Light> transform(const AffineSpace3f& xfm,
                         light_mask_t illumMask,
                         light_mask_t shadowMask) const {
      return new AmbientLight(L,illumMask,shadowMask);
    }

    Color Le(const Vector3f& wo) const {
      return L;
    }

    Color eval(const DifferentialGeometry& dg, const Vector3f& wi) const {
      return L;
    }

    Color sample(const DifferentialGeometry& dg, Sample3f& wi, float& tMax, const Vec2f& s) const {
      wi = cosineSampleHemisphere(s.x, s.y, dg.Ns);
      tMax = inf;
      return L;
    }

    float pdf(const DifferentialGeometry& dg, const Vector3f& wi) const {
      return cosineSampleHemispherePDF(wi,dg.Ns);
    }

  protected:
    Color L;          //!< Radiance (W/(m^2*sr))
  };
}

#endif
