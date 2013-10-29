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

#ifndef __EMBREE_MICROFACET_FRESNEL_H__
#define __EMBREE_MICROFACET_FRESNEL_H__

#include "../../brdfs/optics.h"

namespace embree
{
  /*! Fresnel term for perfect reflection. */
  class FresnelNone {
  public:

    /*! Evaluation simply returns white. */
    __forceinline Color eval(float cosTheta) const { return one; }
  };

  /*! Constant fresnel term for perfect reflection. */
  class FresnelConstant {
  public:

    /*! Constant fresnel term constructor. \param c is the constant
     *  reflectivity to return. */
    __forceinline FresnelConstant(Color c) : c(c) {}

    /*! Evaluation simply returns constant. */
    __forceinline Color eval(float cosTheta) const { return c; }

  private:
    Color c;  //!< constant reflectivity
  };

  /*! Schlick approximation of fresnel term of a dielectric
   *  surface. */
  class SchlickDielectric {
  public:

    /*! Fresnel term constructor. \param f0 is the reflectivity of the
     *  material along the normal vector. */
    __forceinline SchlickDielectric(float f0) : f0(f0) {}

    /*! Evaluates the fresnel term. \param cosTheta is the cosine
     *  between the facet normal (half vector) and the viewing
     *  vector. */
    __forceinline Color eval(float cosTheta) const {
      const float b = 1.0f-cosTheta;
      return Color(f0+(1.0f-f0)*b*b*b*b*b);
    }

  private:
    const float f0;  //!< reflectivity of the material along the normal vector
  };

  /*! Fresnel term of a dielectric surface. A dielectric surface is
   *  for instance glass or water. */
  class FresnelDielectric {
  public:

    /*! Dielectric fresnel term constructor. \param etai is the
     *  refraction index of the medium the incident ray travels in
     *  \param etat is the refraction index of the opposite medium */
    __forceinline FresnelDielectric(float etai, float etat) : etai(etai), etat(etat) {}

    /*! Evaluates the fresnel term. \param cosTheta is the cosine
     *  between the facet normal (half vector) and the viewing
     *  vector. */
    __forceinline Color eval(float cosTheta) const {
      return Color(fresnelDielectric(cosTheta,etai*rcp(etat)));
    }

  private:

    /*! refraction index of the medium the incident ray travels in */
    const float etai;

    /*! refraction index of the medium the outgoing transmission rays
     *  travels in */
    const float etat;
  };

  /*! Fresnel term for a metal surface. */
  class FresnelConductor {
  public:

    /*! Conductor fresnel term constructor. \param eta is the real part of
     *  the refraction index \param k is the imaginary part of the
     *  refraction index */
    __forceinline FresnelConductor(const Color& eta, const Color& k) : eta(eta), k(k) {}

    /*! Evaluates the fresnel term. \param cosTheta is the cosine
     *  between the facet normal (half vector) and the viewing
     *  vector. */
    __forceinline Color eval(float cosTheta) const { return fresnelConductor(cosTheta,eta,k); }
  private:
    const Color eta;  //!< Real part of refraction index
    const Color k;    //!< Imaginary part of refraction index
  };
}

#endif
