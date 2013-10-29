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

#ifndef __EMBREE_BRUSHED_METAL_H__
#define __EMBREE_BRUSHED_METAL_H__

#include "../materials/material.h"
#include "../brdfs/conductor.h"
#include "../brdfs/microfacet.h"

namespace embree
{
  /*! Implements a brushed metal material. The metal is modelled by a
   *  microfacet BRDF with a fresnel term for metals and a power
   *  cosine distribution for different roughness of the material. */
  class BrushedMetal : public Material
  {
    typedef Microfacet<FresnelConductor,AnisotropicPowerCosineDistribution > MicrofacetMetal1;
    typedef Microfacet<FresnelConductor,AnisotropicBeckmannDistribution > MicrofacetMetal2;

  public:

    /*! Construction from parameters. */
    BrushedMetal (const Parms& parms) {
      reflectance  = parms.getColor("reflectance",one);
      eta          = parms.getColor("eta",Color(1.4f));
      k            = parms.getColor("k",Color(0.0f));
      roughnessX   = parms.getFloat("roughnessX",0.01f);
      roughnessY   = parms.getFloat("roughnessY",0.01f);
      rcpRoughnessX = rcp(roughnessX);
      rcpRoughnessY = rcp(roughnessY);
    }

    void shade(const Ray& ray, const Medium& currentMedium, const DifferentialGeometry& dg, CompositedBRDF& brdfs) const
    {
      /*! handle the special case of a specular metal through a special BRDF */
      if (unlikely(roughnessX == 0.0f || roughnessY == 0.0f)) {
        brdfs.add(NEW_BRDF(Conductor)(reflectance, eta, k));
        return;
      }

      brdfs.add(NEW_BRDF(MicrofacetMetal1)(reflectance, FresnelConductor(eta,k), AnisotropicPowerCosineDistribution(dg.Tx,rcpRoughnessX,dg.Ty,rcpRoughnessY,dg.Ns)));
      //brdfs.add(NEW_BRDF(MicrofacetMetal2)(reflectance, FresnelConductor(eta,k), AnisotropicBeckmannDistribution(dx,roughnessX,dy,inf,dz)));
    }

  private:
    Color reflectance;   //!< Reflectivity of the metal
    Color eta;           //!< Real part of refraction index
    Color k;             //!< Imaginary part of refraction index
    float roughnessX;    //!< Roughness parameter in X direction. The range goes from 0 (specular) to 1 (diffuse).
    float roughnessY;    //!< Roughness parameter in Y direction. The range goes from 0 (specular) to 1 (diffuse).
    float rcpRoughnessX; //!< Reciprocal roughness parameter in X direction.
    float rcpRoughnessY; //!< Reciprocal roughness parameter in Y direction.
  };
}

#endif
