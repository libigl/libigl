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

#ifndef __EMBREE_SPOT_LIGHT_H__
#define __EMBREE_SPOT_LIGHT_H__

#include "../lights/light.h"

namespace embree
{
  /*! Implements a spot light source. */
  class SpotLight : public Light
  {
    /*! Construction from position, direction, intensity and opening angles. */
    SpotLight(const Vector3f& P, const Vector3f& _D, const Color& I, 
              float cosAngleMin, float cosAngleMax,
              light_mask_t illumMask=-1,
              light_mask_t shadowMask=-1)
      : Light(illumMask,shadowMask), 
        P(P), _D(_D), I(I), cosAngleMin(cosAngleMin), cosAngleMax(cosAngleMax) 
    {}

  public:

    /*! Construction from parameter container. */
    SpotLight (const Parms& parms) {
      P           = parms.getVector3f("P");
      _D          = -normalize(parms.getVector3f("D"));
      I           = parms.getColor("I");
      cosAngleMin = cosf(0.5f*deg2rad(parms.getFloat("angleMin")));
      cosAngleMax = cosf(0.5f*deg2rad(parms.getFloat("angleMax")));
    }

    Ref<Light> transform(const AffineSpace3f& xfm,
                         light_mask_t illumMask,
                         light_mask_t shadowMask) const {
      return new SpotLight(xfmPoint(xfm,P),xfmVector(xfm,_D),I,
                           cosAngleMin,cosAngleMax,
                           illumMask,shadowMask);
    }

    Color sample(const DifferentialGeometry& dg, Sample3f& wi, float& tMax, const Vec2f& s) const
    {
      Vector3f d = P - dg.P;
      float distance = length(d);
      wi = Sample3f(d * rcp(distance), distance*distance);
      tMax = distance;
      float cosAngle = dot(wi.value,_D);
      if (cosAngleMin != cosAngleMax)
        return I*clamp((cosAngle - cosAngleMax)*rcp(cosAngleMin - cosAngleMax));
      else if (cosAngle > cosAngleMin)
        return I;
      else 
        return zero;
    }

    float pdf (const DifferentialGeometry& dg, const Vector3f& wi) const {
      return zero;
    }

  private:
    Vector3f P;                        //!< Position of the spot light
    Vector3f _D;                       //!< Negative light direction of the spot light
    Color I;                        //!< Radiant intensity (W/sr)
    float cosAngleMin, cosAngleMax; //!< Linear falloff region
  };
}

#endif
