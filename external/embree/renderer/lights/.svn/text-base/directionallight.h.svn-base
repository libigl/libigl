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

#ifndef __EMBREE_DIRECTIONAL_LIGHT_H__
#define __EMBREE_DIRECTIONAL_LIGHT_H__

#include "lights/light.h"

namespace embree
{
  /*! Implements a directional light. */
  class DirectionalLight : public Light
  {
    /*! Construction from members. */
    DirectionalLight (const Vec3f& _wo, const Col3f& E)
      : _wo(normalize(_wo)), E(E) { }

  public:

    /*! Construction from parameter container. */
    DirectionalLight (const Parms& parms) {
      _wo = -normalize(parms.getVec3f("D"));
      E   = parms.getCol3f("E");
    }

    Ref<Light> transform(const AffineSpace3f& xfm) const {
      return new DirectionalLight(xfmVector(xfm,_wo),E);
    }

    Col3f sample(const DifferentialGeometry& dg, Sample3f& wi, float& tMax, const Vec2f& s) const {
      wi = _wo; tMax = inf; return E;
    }

    float pdf(const DifferentialGeometry& dg, const Vec3f& wi) const {
      return zero;
    }

  private:
    Vec3f _wo;    //!< inverse light direction
    Col3f E;      //!< Irradiance (W/m^2)
  };
}

#endif
