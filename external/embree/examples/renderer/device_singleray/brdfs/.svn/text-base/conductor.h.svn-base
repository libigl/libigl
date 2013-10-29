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

#ifndef __EMBREE_CONDUCTOR_BRDF_H__
#define __EMBREE_CONDUCTOR_BRDF_H__

#include "../brdfs/brdf.h"
#include "../brdfs/optics.h"

namespace embree
{
  /*! BRDF of a perfectly polished metal. */
  class Conductor : public BRDF
  {
  public:

    /*! Conductor BRDF constructor. \param R is the reflectivity
     *  coefficient \param eta the real part of the refraction index
     *  and \param k the imaginary part of the refraction index */
    __forceinline Conductor(const Color& R, const Color& eta, const Color& k)
      : BRDF(SPECULAR_REFLECTION), R(R), eta(eta), k(k) {}

    __forceinline Color eval(const Vector3f& wo, const DifferentialGeometry& dg, const Vector3f& wi) const {
      return zero;
    }

    Color sample(const Vector3f& wo, const DifferentialGeometry& dg, Sample3f& wi, const Vec2f& s) const {
      wi = reflect(wo,dg.Ns);
      return R * fresnelConductor(dot(wo,dg.Ns),eta,k);
    }

    float pdf(const Vector3f& wo, const DifferentialGeometry& dg, const Vector3f& wi) const {
      return zero;
    }

  private:
    Color R;          //!< Reflectivity coefficient
    Color eta;        //!< Real part of refraction index
    Color k;          //!< Imaginary part of refraction index
  };
}

#endif
