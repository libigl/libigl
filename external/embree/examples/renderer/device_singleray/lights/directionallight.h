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

#ifndef __EMBREE_DIRECTIONAL_LIGHT_H__
#define __EMBREE_DIRECTIONAL_LIGHT_H__

#include "../lights/light.h"

namespace embree
{
  /*! Implements a directional light. */
  class DirectionalLight : public Light
  {
    /*! Construction from members. */
    DirectionalLight (const Vector3f& _wo, const Color& E,
                      light_mask_t illumMask=-1,
                      light_mask_t shadowMask=-1)
      : Light(illumMask,shadowMask), 
        _wo(normalize(_wo)), E(E) { }

  public:

    /*! Construction from parameter container. */
    DirectionalLight (const Parms& parms) {
      _wo = -normalize(parms.getVector3f("D"));
      E   = parms.getColor("E");
    }

    Ref<Light> transform(const AffineSpace3f& xfm,
                         light_mask_t illumMask,
                         light_mask_t shadowMask) const {
      return new DirectionalLight(xfmVector(xfm,_wo),E,illumMask,shadowMask);
    }

    Color sample(const DifferentialGeometry& dg, Sample3f& wi, float& tMax, const Vec2f& s) const {
      wi = _wo; tMax = inf; return E;
    }

    float pdf(const DifferentialGeometry& dg, const Vector3f& wi) const {
      return zero;
    }

  private:
    Vector3f _wo;    //!< negative light direction
    Color E;      //!< Irradiance (W/m^2)
  };
}

#endif
