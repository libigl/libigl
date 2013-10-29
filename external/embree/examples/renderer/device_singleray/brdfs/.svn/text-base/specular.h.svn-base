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

#ifndef __EMBREE_SPECULAR_BRDF_H__
#define __EMBREE_SPECULAR_BRDF_H__

#include "../brdfs/brdf.h"
#include "../brdfs/optics.h"

namespace embree
{
  /*! Specular Phong BRDF. The cosine of the angle between incoming
   *  ray direction and reflection direction is raised to some power
   *  to approximate a glossy reflection. */
  class Specular : public BRDF
  {
  public:

    /*! Specular BRDF constructor. */
    __forceinline Specular(const Color& R, float exp) : BRDF(GLOSSY_REFLECTION), R(R), exp(exp) {}

    __forceinline Color eval(const Vector3f& wo, const DifferentialGeometry& dg, const Vector3f& wi) const {
      Vector3f r = reflect(wo,dg.Ns);
      if (dot(r,wi) < 0) return zero;
      return R * (exp+2) * (1.0f/(2.0f*float(pi))) * pow(dot(r,wi),exp) * clamp(dot(wi,dg.Ns));
    }

    Color sample(const Vector3f& wo, const DifferentialGeometry& dg, Sample3f& wi, const Vec2f& s) const {
      return eval(wo, dg, wi = powerCosineSampleHemisphere(s.x,s.y,reflect(wo,dg.Ns),exp));
    }

    float pdf(const Vector3f& wo, const DifferentialGeometry& dg, const Vector3f& wi) const {
      return powerCosineSampleHemispherePDF(wi,reflect(wo,dg.Ns),exp);
    }

  private:

    /*! The reflectivity parameter. The range is [0,1] where 0 means
     *  no reflection at all, and 1 means full reflection. */
    Color R;

    /*! The exponent that determines the glossiness. The range is
     *  [0,infinity[ where 0 means a diffuse surface, and the
     *  specularity increases towards infinity. */
    float exp;
  };
}

#endif
