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

#ifndef __EMBREE_DISTANT_LIGHT_H__
#define __EMBREE_DISTANT_LIGHT_H__

#include "../lights/light.h"

namespace embree
{
  /*! Implements a distant light. The distant light illuminates from
   *  infinity from a cone of directions. This simulates the light
   *  field of a far away big object like the sun. */
  class DistantLight : public EnvironmentLight
  {
  protected:

    /*! Construction from negative light direction, radiance and
     *  opening angle of cone. */
    DistantLight (const Vector3f& _wo, const Color& L, float halfAngle,
                  light_mask_t illumMask=-1,
                  light_mask_t shadowMask=-1)
      : EnvironmentLight(illumMask,shadowMask), 
        _wo(normalize(_wo)), L(L), halfAngle(halfAngle), 
        cosHalfAngle(cosf(halfAngle)) 
    {}

  public:

    /*! Construction from parameter container. */
    DistantLight (const Parms& parms) {
      _wo       = -normalize(parms.getVector3f("D"));
      L         = parms.getColor("L");
      halfAngle = deg2rad(parms.getFloat("halfAngle"));
      cosHalfAngle = cosf(halfAngle);
    }

    Ref<Light> transform(const AffineSpace3f& xfm,
                         light_mask_t illumMask,
                         light_mask_t shadowMask) const {
      return new DistantLight(xfmVector(xfm,_wo),L,halfAngle,
                              illumMask,shadowMask);
    }

    Color Le(const Vector3f& wo) const {
      if (dot(-wo,_wo) >= cosHalfAngle) return L;
      return zero;
    }

    Color eval(const DifferentialGeometry& dg, const Vector3f& wi) const {
      if (dot(wi,_wo) >= cosHalfAngle) return L;
      return zero;
    }

    Color sample(const DifferentialGeometry& dg, Sample3f& wi, float& tMax, const Vec2f& s) const {
      wi = uniformSampleCone(s.x,s.y,halfAngle,_wo);
      tMax = inf;
      return L;
    }

    float pdf(const DifferentialGeometry& dg, const Vector3f& wi) const {
      return uniformSampleConePDF(wi,halfAngle,_wo);
    }

  protected:
    Vector3f _wo;           //!< Negative light direction
    Color L;             //!< Radiance (W/(m^2*sr))
    float halfAngle;     //!< Half illumination angle
    float cosHalfAngle;  //!< Cosine of half illumination angle
  };
}

#endif
