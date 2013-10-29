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

#ifndef __EMBREE_MICROFACET_BECKMANN_DISTRIBUTION_H__
#define __EMBREE_MICROFACET_BECKMANN_DISTRIBUTION_H__

#include "../../default.h"

namespace embree
{
  /*! Beckmann microfacet distribution. */
  class BeckmannDistribution {
  public:

    /*! Distribution constructor. */
    __forceinline BeckmannDistribution(float sigma, const Vector3f& dz) 
      : sigma(sigma), rcpSigma(rcp(sigma)), 
        dz(dz), norm(float(pi)*sigma*sigma) {}

    /*! Evaluates the distribution. \param wh is the half vector */
    __forceinline float eval(const Vector3f& wh) const 
    {
      const float cosTheta = dot(wh,dz);
      const float cosTheta2 = sqr(cosTheta);
      const float sinTheta2 = 1-cosTheta2;
      return exp(-sqr(rcpSigma)*sinTheta2*rcp(cosTheta2))*rcp(norm*sqr(cosTheta2));
    }

    /*! Samples the distribution. \param s is the sample location
     *  provided by the caller. */
    __forceinline Sample3f sample(const Vec2f& s) const
    {
      const float phi = float(two_pi)*s.x;
      const float cosPhi = cosf(phi);
      const float sinPhi = sinf(phi);
      const float rcpSigma2 = sqr(rcpSigma);
      const float sinThetaT = sqrt(-log(s.y));
      const float cosThetaT = sqrt(rcpSigma2);
      const float rcpT = rsqrt(sqr(cosThetaT)+sqr(sinThetaT));
      const float cosTheta = cosThetaT*rcpT, cosTheta2 = sqr(cosTheta);
      const float sinTheta = sinThetaT*rcpT, sinTheta2 = sqr(sinTheta);
      const float pdf = exp(-rcpSigma2*sinTheta2*rcp(cosTheta2))*rcp(norm*cosTheta*cosTheta2);
      const Vector3f wh = frame(dz)*Vector3f(cosPhi*sinTheta,sinPhi*sinTheta,cosTheta);
      return Sample3f(wh,pdf);
    }

    /*! Evaluates the sampling PDF. \param wh is the direction to
     *  evaluate the PDF for \returns the probability density */
    __forceinline float pdf(const Vector3f& wh)  const
    {
      const float cosTheta = dot(wh,dz);
      if (cosTheta < 0.0f) return 0.0f;
      const float cosTheta2 = sqr(cosTheta);
      const float sinTheta2 = 1-cosTheta2;
      return exp(-sqr(rcpSigma)*sinTheta2*rcp(cosTheta2))*rcp(norm*cosTheta*cosTheta2);
    }

  private:
    const float sigma;      //!< Standard deviation of distribution.
    const float rcpSigma;   //!< Reciprocal standard deviation of distribution.
    const Vector3f dz;         //!< z-direction of the distribution.
    const float norm;       //!< Normalization constant.
  };
}

#endif
