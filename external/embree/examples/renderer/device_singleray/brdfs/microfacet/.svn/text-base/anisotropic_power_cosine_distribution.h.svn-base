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

#ifndef __EMBREE_MICROFACET_ANISOTROPIC_POWER_COSINE_DISTRIBUTION_H__
#define __EMBREE_MICROFACET_ANISOTROPIC_POWER_COSINE_DISTRIBUTION_H__

#include "../../default.h"

namespace embree
{
  /*! Anisotropic power cosine microfacet distribution. */
  class AnisotropicPowerCosineDistribution {
  public:

    /*! Anisotropic power cosine distribution constructor. */
    __forceinline AnisotropicPowerCosineDistribution(const Vector3f& dx, float nx, const Vector3f& dy, float ny, const Vector3f& dz) 
      : dx(dx), nx(nx), dy(dy), ny(ny), dz(dz),
        norm1(sqrtf((nx+1)*(ny+1)) * float(one_over_two_pi)),
        norm2(sqrtf((nx+2)*(ny+2)) * float(one_over_two_pi)) {}

    /*! Evaluates the power cosine distribution. \param wh is the half
     *  vector */
    __forceinline float eval(const Vector3f& wh) const 
    {
      const float cosPhiH   = dot(wh, dx);
      const float sinPhiH   = dot(wh, dy);
      const float cosThetaH = dot(wh, dz);
      const float R = sqr(cosPhiH)+sqr(sinPhiH);
      if (R == 0.0f) return norm2;
      const float n = (nx*sqr(cosPhiH)+ny*sqr(sinPhiH))*rcp(R);
      return norm2 * pow(abs(cosThetaH), n);
    }

    /*! Samples the distribution. \param s is the sample location
     *  provided by the caller. */
    __forceinline Sample3f sample(const Vec2f& s) const
    {
      const float phi = float(two_pi)*s.x;
      const float sinPhi0 = sqrtf(nx+1)*sinf(phi);
      const float cosPhi0 = sqrtf(ny+1)*cosf(phi);
      const float norm = rsqrt(sqr(sinPhi0)+sqr(cosPhi0));
      const float sinPhi = sinPhi0*norm;
      const float cosPhi = cosPhi0*norm;
      const float n = nx*sqr(cosPhi)+ny*sqr(sinPhi);
      const float cosTheta = powf(s.y,rcp(n+1));
      const float sinTheta = cos2sin(cosTheta);
      const float pdf = norm1*powf(cosTheta,n);
      const Vector3f h(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta);
      const Vector3f wh = h.x*dx + h.y*dy + h.z*dz;
      return Sample3f(wh,pdf);
    }

    /*! Evaluates the sampling PDF. \param wh is the direction to
     *  evaluate the PDF for \returns the probability density */
    __forceinline float pdf(const Vector3f& wh) const
                            
    {
      const float cosPhiH   = dot(wh, dx);
      const float sinPhiH   = dot(wh, dy);
      const float cosThetaH = dot(wh, dz);
      const float R = sqr(cosPhiH)+sqr(sinPhiH);
      if (R == 0.0f) return norm1;
      const float n = (nx*sqr(cosPhiH)+ny*sqr(sinPhiH))*rcp(R);
      return norm1*powf(cosThetaH,n);
    }

  private:
    const Vector3f dx;     //!< x-direction of the distribution.
    const float nx;     //!< Glossiness in x direction with range [0,infinity[ where 0 is a diffuse surface.
    const Vector3f dy;     //!< y-direction of the distribution.
    const float ny;     //!< Exponent that determines the glossiness in y direction.
    const Vector3f dz;     //!< z-direction of the distribution.
    const float norm1;  //!< Normalization constant for calculating the pdf for sampling.
    const float norm2;  //!< Normalization constant for calculating the distribution.
  };
}

#endif
