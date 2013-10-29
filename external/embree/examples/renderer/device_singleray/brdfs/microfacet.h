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

#ifndef __EMBREE_MICROFACET_BRDF_H__
#define __EMBREE_MICROFACET_BRDF_H__

#include "../brdfs/brdf.h"
#include "../brdfs/optics.h"
#include "microfacet/fresnel.h"

#include "../brdfs/microfacet/power_cosine_distribution.h"
#include "../brdfs/microfacet/anisotropic_power_cosine_distribution.h"

#include "../brdfs/microfacet/beckmann_distribution.h"
#include "../brdfs/microfacet/anisotropic_beckmann_distribution.h"

namespace embree
{
  /*! Microfacet BRDF model. The class is templated over the fresnel
   *  term and microfacet distribution to use. */
  template<typename Fresnel, typename Distribution>
    class Microfacet : public BRDF
  {
  public:

    /*! Microfacet BRDF constructor. This is a glossy BRDF. */
    __forceinline Microfacet (const Color& R, const Fresnel& fresnel, const Distribution& distribution)
    : BRDF(GLOSSY_REFLECTION), R(R), fresnel(fresnel), distribution(distribution) {}
    
    __forceinline Color eval(const Vector3f& wo, const DifferentialGeometry& dg, const Vector3f& wi) const
    {
      if (dot(wi,dg.Ng) <= 0) return zero;
      const float cosThetaO = dot(wo,dg.Ns);
      const float cosThetaI = dot(wi,dg.Ns);
      if (cosThetaI <= 0.0f || cosThetaO <= 0.0f) return zero;
      const Vector3f wh = normalize(wi + wo);
      const float cosThetaH = dot(wh, dg.Ns);
      const float cosTheta = dot(wi, wh); // = dot(wo, wh);
      const Color F = fresnel.eval(cosTheta);
      const float D = distribution.eval(wh);
      const float G = min(1.0f, 2.0f * cosThetaH * cosThetaO * rcp(cosTheta), 2.0f * cosThetaH * cosThetaI * rcp(cosTheta));
      return R * D * G * F * rcp(4.0f*cosThetaO);
    }

    Color sample(const Vector3f& wo, const DifferentialGeometry& dg, Sample3f& wi, const Vec2f& s) const
    {
      if (dot(wo,dg.Ns) <= 0.0f) return zero;
      const Sample3f wh = distribution.sample(s);
      wi = Sample3f(reflect(wo,wh),wh.pdf*rcp(4.0f*abs(dot(wo,(Vector3f)wh))));
      if (dot((Vector3f)wi,dg.Ns) <= 0.0f) return zero;
      return eval(wo,dg,wi);
    }

    float pdf(const Vector3f& wo, const DifferentialGeometry& dg, const Vector3f& wi) const {
      const Vector3f wh = normalize(wo+wi);
      return distribution.pdf(wh)*rcp(4.0f*abs(dot(wo,wh)));
    }

  private:

    /*! Reflectivity of the microfacets. The range is [0,1] where 0
     *  means no reflection at all, and 1 means full reflection. */
    const Color R;

    /*! Fresnel term to use. */
    const Fresnel fresnel;

    /*! Microfacet distribution to use. */
    const Distribution distribution;
  };
}

#endif
