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

#ifndef __EMBREE_METAL_H__
#define __EMBREE_METAL_H__

#include "../materials/material.h"
#include "../brdfs/conductor.h"
#include "../brdfs/microfacet.h"

namespace embree
{
  /*! Implements a rough metal material. The metal is modelled by a
   *  microfacet BRDF with a fresnel term for metals and a power
   *  cosine distribution for different roughness of the material. */
  class Metal : public Material
  {
    typedef Microfacet<FresnelConductor,PowerCosineDistribution > MicrofacetMetal;
  public:

    /*! Construction from parameters. */
    Metal (const Parms& parms) {
      reflectance  = parms.getColor("reflectance",one);
      eta          = parms.getColor("eta",Color(1.4f));
      k            = parms.getColor("k",Color(0.0f));
      roughness    = parms.getFloat("roughness",0.01f);
      rcpRoughness = rcp(roughness);
    }

    void shade(const Ray& ray, const Medium& currentMedium, const DifferentialGeometry& dg, CompositedBRDF& brdfs) const
    {
      /*! handle the special case of a specular metal through a special BRDF */
      if (roughness == 0.0f)
        brdfs.add(NEW_BRDF(Conductor)(reflectance, eta, k));

      /*! otherwise use the microfacet BRDF model */
      else
        brdfs.add(NEW_BRDF(MicrofacetMetal)(reflectance, FresnelConductor(eta,k), PowerCosineDistribution(rcpRoughness,dg.Ns)));
    }

  protected:
    Color reflectance; //!< Reflectivity of the metal
    Color eta;         //!< Real part of refraction index
    Color k;           //!< Imaginary part of refraction index
    float roughness;   //!< Roughness parameter. The range goes from 0 (specular) to 1 (diffuse).
    float rcpRoughness;//!< Reciprocal roughness parameter.
  };
}

#endif
