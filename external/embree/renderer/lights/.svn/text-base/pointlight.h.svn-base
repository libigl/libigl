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

#ifndef __EMBREE_POINT_LIGHT_H__
#define __EMBREE_POINT_LIGHT_H__

#include "lights/light.h"

namespace embree
{
  /*! Implements a point light source. */
  class PointLight : public Light
  {
    /*! Construction from position and intensity. */
    PointLight (const Vec3f& P, const Col3f& I) : P(P), I(I) {}

  public:

    /*! Construction from parameter container. */
    PointLight (const Parms& parms) {
      P = parms.getVec3f("P",zero);
      I = parms.getCol3f("I",zero);
    }

    Ref<Light> transform(const AffineSpace3f& xfm) const {
      return new PointLight(xfmPoint(xfm,P),I);
    }

    Col3f sample(const DifferentialGeometry& dg, Sample3f& wi, float& tMax, const Vec2f& s) const
    {
      Vec3f d = P - dg.P;
      float distance = length(d);
      wi = Sample3f(d / distance, distance*distance);
      tMax = distance;
      return I;
    }

    float pdf(const DifferentialGeometry& dg, const Vec3f& wi) const {
      return zero;
    }

  private:
    Vec3f P;       //!< Position of the point light
    Col3f I;       //!< Radiant intensity (W/sr)
  };
}

#endif
