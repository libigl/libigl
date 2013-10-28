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

#ifndef __EMBREE_MICROFACET_ANISOTROPIC_BECKMANN_DISTRIBUTION_H__
#define __EMBREE_MICROFACET_ANISOTROPIC_BECKMANN_DISTRIBUTION_H__

#include "../../default.h"

namespace embree
{
  /*! Anisotropic Beckmann microfacet distribution. */
  class AnisotropicBeckmannDistribution {
  public:

    /*! Distribution constructor. */
    __forceinline AnisotropicBeckmannDistribution(const Vector3f& dx, float sigmaX, const Vector3f& dy, float sigmaY, const Vector3f& dz) 
      : dx(dx), sigmaX(sigmaX), rcpSigmaX(rcp(sigmaX)), 
        dy(dy), sigmaY(sigmaY), rcpSigmaY(rcp(sigmaY)),
        dz(dz), norm(float(pi)*sigmaX*sigmaY) {}

    /*! Evaluates the distribution. \param wh is the half vector */
    __forceinline float eval(const Vector3f& wh) const 
    {
      const float cosPhiR = dot(wh,dx);
      const float sinPhiR = dot(wh,dy);
      const float R = sqr(cosPhiR)+sqr(sinPhiR);
      const float cosTheta = dot(wh,dz);
      const float cosTheta2 = sqr(cosTheta);
      const float sinTheta2 = 1-cosTheta2;
      if (R == 0.0f) return rcp(norm*sqr(cosTheta2));
      const float rcpSigma2 = (sqr(cosPhiR*rcpSigmaX)+sqr(sinPhiR*rcpSigmaY))*rcp(R);
      return exp(-rcpSigma2*sinTheta2*rcp(cosTheta2))*rcp(norm*sqr(cosTheta2));
    }

    /*! Samples the distribution. \param s is the sample location
     *  provided by the caller. */
    __forceinline Sample3f sample(const Vec2f& s) const
    {
      const float phi = float(two_pi)*s.x;
      const float sinPhiR = sigmaY*sinf(phi);
      const float cosPhiR = sigmaX*cosf(phi);
      const float rcpR = rsqrt(sqr(sinPhiR)+sqr(cosPhiR));
      const float sinPhi = sinPhiR*rcpR;
      const float cosPhi = cosPhiR*rcpR;
      const float rcpSigma2 = sqr(cosPhi*rcpSigmaX)+sqr(sinPhi*rcpSigmaY);
      const float sinThetaT = sqrt(-log(s.y));
      const float cosThetaT = sqrt(rcpSigma2);
      const float rcpT = rsqrt(sqr(cosThetaT)+sqr(sinThetaT));
      const float cosTheta = cosThetaT*rcpT, cosTheta2 = sqr(cosTheta);
      const float sinTheta = sinThetaT*rcpT, sinTheta2 = sqr(sinTheta);
      const float pdf = exp(-rcpSigma2*sinTheta2*rcp(cosTheta2))*rcp(norm*cosTheta*cosTheta2);
      const Vector3f wh = cosPhi*sinTheta*dx + sinPhi*sinTheta*dy + cosTheta*dz;
      return Sample3f(wh,pdf);
    }

    /*! Evaluates the sampling PDF. \param wh is the direction to
     *  evaluate the PDF for \returns the probability density */
    __forceinline float pdf(const Vector3f& wh)  const
    {
      const float cosPhiR  = dot(wh, dx);
      const float sinPhiR  = dot(wh, dy);
      const float cosTheta = dot(wh, dz);
      if (cosTheta < 0.0f) return 0.0f;
      const float R = sqr(cosPhiR)+sqr(sinPhiR);
      if (R == 0.0f) return rcp(norm*cosTheta*sqr(cosTheta));
      const float rcpSigma2 = (sqr(cosPhiR*rcpSigmaX)+sqr(sinPhiR*rcpSigmaY))*rcp(R);
      const float cosTheta2 = sqr(cosTheta);
      const float sinTheta2 = 1-cosTheta2;
      return exp(-rcpSigma2*sinTheta2*rcp(cosTheta2))*rcp(norm*cosTheta*cosTheta2);
    }

  private:
    const Vector3f dx;         //!< x-direction of the distribution.
    const float sigmaX;     //!< Standard deviation of distribution in X direction.
    const float rcpSigmaX;  //!< Reciprocal standard deviation of distribution in X direction.
    const Vector3f dy;         //!< y-direction of the distribution.
    const float sigmaY;     //!< Standard deviation of distribution in Y direction.
    const float rcpSigmaY;  //!< Reciprocal standard deviation of distribution in X direction.
    const Vector3f dz;         //!< z-direction of the distribution.
    const float norm;       //!< Normalization constant.
  };
}

#endif
