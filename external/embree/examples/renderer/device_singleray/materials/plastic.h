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

#ifndef __EMBREE_PLASTIC_H__
#define __EMBREE_PLASTIC_H__

#include "../materials/material.h"
#include "../brdfs/lambertian.h"
#include "../brdfs/dielectric.h"
#include "../brdfs/dielectriclayer.h"
#include "../brdfs/microfacet.h"

namespace embree
{
  /*! Implements a plastic material. A dielectric layer over a diffuse
   *  surface is modeled through a dieletric layer BRDF for the
   *  covered diffuse part, and a microfacet BRDF for the rough
   *  dielectric surface. */
  class Plastic : public Material
  {
    typedef Microfacet<FresnelDielectric,PowerCosineDistribution > MicrofacetPlastic;

  public:

    /*! Construction from parameters. */
    Plastic (const Parms& parms) {
      pigmentColor = parms.getColor("pigmentColor",one);
      eta          = parms.getFloat("eta",1.4f);
      roughness    = parms.getFloat("roughness",0.01f);
      rcpRoughness = rcp(roughness);
    }

    void shade(const Ray& ray, const Medium& currentMedium, const DifferentialGeometry& dg, CompositedBRDF& brdfs) const
    {
      /*! the dielectric layer that models the covered diffuse part */
      brdfs.add(NEW_BRDF(DielectricLayer<Lambertian >)(one, 1.0f, eta, Lambertian (pigmentColor)));

      /*! use dielectric reflection in case of a specular surface */
      if (roughness == 0.0f)
        brdfs.add(NEW_BRDF(DielectricReflection)(1.0f, eta));

      /*! otherwise use the microfacet BRDF to model the rough surface */
      else
        brdfs.add(NEW_BRDF(MicrofacetPlastic)(one, FresnelDielectric(1.0f, eta), PowerCosineDistribution(rcpRoughness,dg.Ns)));
    }

  protected:
    Color pigmentColor; //!< Color of the diffuse layer.
    float eta;          //!< Refraction index of the dielectric layer.
    float roughness;    //!< Roughness parameter. The range goes from 0 (specular) to 1 (diffuse).
    float rcpRoughness; //!< Reciprocal roughness parameter.
  };
}

#endif
