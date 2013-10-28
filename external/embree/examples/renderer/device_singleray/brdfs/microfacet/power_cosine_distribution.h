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

#ifndef __EMBREE_MICROFACET_POWER_COSINE_DISTRIBUTION_H__
#define __EMBREE_MICROFACET_POWER_COSINE_DISTRIBUTION_H__

#include "../../default.h"

namespace embree
{
  /*! Power cosine microfacet distribution. */
  class PowerCosineDistribution {
  public:

    /*! Power cosine distribution constructor. */
    __forceinline PowerCosineDistribution(float n, const Vector3f& dz) 
      : n(n), dz(dz), norm1((n+1)*float(one_over_two_pi)), norm2((n+2)*float(one_over_two_pi)) {}

    /*! Evaluates the power cosine distribution. \param wh is the half
     *  vector */
    __forceinline float eval(const Vector3f& wh) const 
    {
      const float cosTheta = dot(wh,dz);
      return norm2 * pow(abs(cosTheta),n);
    }

    /*! Samples the distribution. \param s is the sample location
     *  provided by the caller. */
    __forceinline Sample3f sample(const Vec2f& s) const
    {
      const float phi = float(two_pi) * s.x;
      const float cosPhi = cos(phi);
      const float sinPhi = sin(phi);
      const float cosTheta = pow(s.y,rcp(n+1));
      const float sinTheta = cos2sin(cosTheta);
      return Sample3f(frame(dz)*Vector3f(cosPhi*sinTheta,sinPhi*sinTheta,cosTheta), norm1*pow(cosTheta,n));
    }

    /*! Evaluates the sampling PDF. \param wh is the direction to
     *  evaluate the PDF for \returns the probability density */
    __forceinline float pdf(const Vector3f& wh) const
    {
      const float cosTheta = dot(wh,dz);
      return norm1*pow(abs(cosTheta),n);
    }

  private:
    const float n;      //!< Glossiness with range [0,infinity[ where 0 is a diffuse surface.
    const Vector3f dz;     //!< z-direction of the distribution.
    const float norm1;  //!< Normalization constant for calculating the pdf for sampling.
    const float norm2;  //!< Normalization constant for calculating the distribution.
  };
}

#endif
