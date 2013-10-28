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

#ifndef __EMBREE_VELVET_H__
#define __EMBREE_VELVET_H__

#include "../materials/material.h"
#include "../brdfs/minnaert.h"
#include "../brdfs/velvety.h"
#include "../brdfs/lambertian.h" 

namespace embree
{
  /*! Implements a velvet material. */
  class Velvet : public Material
  {
  public:

    /*! Construction from parameters. */
    Velvet (const Parms& parms) {
      reflectance = parms.getColor("reflectance",one);
      backScattering = parms.getFloat("backScattering",zero);
      horizonScatteringColor = parms.getColor("horizonScatteringColor",one);
      horizonScatteringFallOff = parms.getFloat("horizonScatteringFallOff",zero);
    }

    void shade(const Ray& ray, const Medium& currentMedium, const DifferentialGeometry& dg, CompositedBRDF& brdfs) const {
      brdfs.add(NEW_BRDF(Minnaert)(reflectance, backScattering));
      brdfs.add(NEW_BRDF(Velvety)(horizonScatteringColor, horizonScatteringFallOff));
    }

  protected:

    /*! Diffuse reflectance of the surface. The range is from 0
     *  (black) to 1 (white). */
    Color reflectance;

    /*! Amount of back scattering. The range is from 0 (no back
     *  scattering) to inf (maximum back scattering). */
    float backScattering;

    /*! Color of horizon scattering. */
    Color horizonScatteringColor;

    /*! Fall-off of horizon scattering. */
    float horizonScatteringFallOff;
  };
}

#endif
