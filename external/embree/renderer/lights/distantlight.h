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

#ifndef __EMBREE_DISTANT_LIGHT_H__
#define __EMBREE_DISTANT_LIGHT_H__

#include "lights/light.h"

namespace embree
{
  /*! Implements a distant light. The distant light illuminates from
   *  infinity from a cone of directions. This simulates the light
   *  field of a far away big object like the sun. */
  class DistantLight : public EnvironmentLight
  {
    /*! Construction from negative light direction, radiance and
     *  opening angle of cone. */
    DistantLight (const Vec3f& _wo, const Col3f& L, float halfAngle)
      : _wo(normalize(_wo)), L(L), halfAngle(halfAngle), cosHalfAngle(cosf(halfAngle)) { }

  public:

    /*! Construction from parameter container. */
    DistantLight (const Parms& parms) {
      _wo       = -normalize(parms.getVec3f("D"));
      L         = parms.getCol3f("L");
      halfAngle = deg2rad(parms.getFloat("halfAngle"));
      cosHalfAngle = cosf(halfAngle);
    }

    Ref<Light> transform(const AffineSpace3f& xfm) const {
      return new DistantLight(xfmVector(xfm,_wo),L,halfAngle);
    }

    Col3f Le(const Vec3f& wo) const {
      if (dot(-wo,_wo) >= cosHalfAngle) return L;
      return zero;
    }

    Col3f eval(const DifferentialGeometry& dg, const Vec3f& wi) const {
      if (dot(wi,_wo) >= cosHalfAngle) return L;
      return zero;
    }

    Col3f sample(const DifferentialGeometry& dg, Sample3f& wi, float& tMax, const Vec2f& s) const {
      wi = uniformSampleCone(s.x,s.y,halfAngle,_wo);
      tMax = inf;
      return L;
    }

    float pdf(const DifferentialGeometry& dg, const Vec3f& wi) const {
      return uniformSampleConePDF(wi,halfAngle,_wo);
    }

  private:
    Vec3f _wo;           //!< Negative light direction
    Col3f L;             //!< Radiance (W/(m^2*sr))
    float halfAngle;     //!< Half illumination angle
    float cosHalfAngle;  //!< Cosine of half illumination angle
  };
}

#endif
