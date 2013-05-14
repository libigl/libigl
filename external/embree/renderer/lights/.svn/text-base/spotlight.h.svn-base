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

#ifndef __EMBREE_SPOT_LIGHT_H__
#define __EMBREE_SPOT_LIGHT_H__

#include "lights/light.h"

namespace embree
{
  /*! Implements a spot light source. */
  class SpotLight : public Light
  {
    /*! Construction from position, direction, intensity and opening angles. */
    SpotLight(const Vec3f& P, const Vec3f& _D, const Col3f& I, float cosAngleMin, float cosAngleMax)
      : P(P), _D(_D), I(I), cosAngleMin(cosAngleMin), cosAngleMax(cosAngleMax) {}

  public:

    /*! Construction from parameter container. */
    SpotLight (const Parms& parms) {
      P           = parms.getVec3f("P");
      _D          = -normalize(parms.getVec3f("D"));
      I           = parms.getCol3f("I");
      cosAngleMin = cosf(deg2rad(parms.getFloat("angleMin")));
      cosAngleMax = cosf(deg2rad(parms.getFloat("angleMax")));
    }

    Ref<Light> transform(const AffineSpace3f& xfm) const {
      return new SpotLight(xfmPoint(xfm,P),xfmVector(xfm,_D),I,cosAngleMin,cosAngleMax);
    }

    Col3f sample(const DifferentialGeometry& dg, Sample3f& wi, float& tMax, const Vec2f& s) const
    {
      Vec3f d = P - dg.P;
      float distance = length(d);
      wi = Sample3f(d * rcp(distance), distance*distance);
      tMax = distance;
      return I*clamp((dot(wi.value,_D) - cosAngleMax)*rcp(cosAngleMin - cosAngleMax));
    }

    float pdf (const DifferentialGeometry& dg, const Vec3f& wi) const {
      return zero;
    }

  private:
    Vec3f P;                        //!< Position of the spot light
    Vec3f _D;                       //!< Negative light direction of the spot light
    Col3f I;                        //!< Radiant intensity (W/sr)
    float cosAngleMin, cosAngleMax; //!< Linear falloff region
  };
}

#endif
