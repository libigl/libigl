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

#ifndef __EMBREE_THIN_DIELECTRIC_H__
#define __EMBREE_THIN_DIELECTRIC_H__

#include "../brdfs/dielectric.h"
#include "../materials/material.h"

namespace embree
{
  /*! Implements a thin dielectricum material. The model uses a
   *  dielectric reflection BRDF and thin dielectric transmission
   *  BRDF. */
  class ThinDielectric : public Material
  {
  public:

    /*! Construction from parameters. */
    ThinDielectric(const Parms& parms) {
      transmission = parms.getColor("transmission",one);
      eta          = parms.getFloat("eta",1.4f);
      thickness    = parms.getFloat("thickness",0.1f);
    }

    void shade(const Ray& ray, const Medium& currentMedium, const DifferentialGeometry& dg, CompositedBRDF& brdfs) const {
      brdfs.add(NEW_BRDF(DielectricReflection)(1.0f, eta));
      brdfs.add(NEW_BRDF(ThinDielectricTransmission)(1.0f, eta, transmission, thickness));
    }

  protected:
    Color transmission;   //!< Transmission coefficient of material.
    float eta;            //!< Refraction index of material.
    float thickness;      //!< Thickness of material layer.
  };
}

#endif
