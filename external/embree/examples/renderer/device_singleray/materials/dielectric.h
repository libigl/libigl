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

#ifndef __EMBREE_DIELECTRIC_H__
#define __EMBREE_DIELECTRIC_H__

#include "../materials/material.h"
#include "../brdfs/dielectric.h"

namespace embree
{
  /*! Implements a dielectric material such as glass. */
  class Dielectric : public Material
  {
  public:

    /*! Construction from parameters. */
    Dielectric(const Parms& parms)
    {
      /*! extract parameters */
      mediumOutside.eta          = parms.getFloat("etaOutside",1.0f);
      mediumInside.eta           = parms.getFloat("etaInside",1.4f);
      mediumOutside.transmission = parms.getColor("transmissionOutside",one);
      mediumInside.transmission  = parms.getColor("transmission",one);
      isMediaInterface           = true;

      /*! precompute BRDF components for more efficient shading */
      reflection_io   = new DielectricReflection  (mediumInside.eta, mediumOutside.eta);
      reflection_oi   = new DielectricReflection  (mediumOutside.eta,mediumInside.eta);
      transmission_io = new DielectricTransmission(mediumInside.eta, mediumOutside.eta);
      transmission_oi = new DielectricTransmission(mediumOutside.eta,mediumInside.eta);
    }

    /*! Destruction */
    ~Dielectric()
    {
      if (reflection_io  ) delete reflection_io;   reflection_io = NULL;
      if (reflection_oi  ) delete reflection_oi;   reflection_oi = NULL;
      if (transmission_io) delete transmission_io; transmission_io = NULL;
      if (transmission_oi) delete transmission_oi; transmission_oi = NULL;
    }

    void shade(const Ray& ray, const Medium& currentMedium, const DifferentialGeometry& dg, CompositedBRDF& brdfs) const
    {
      /*! the ray transitions from outside to inside */
      if (currentMedium == mediumOutside) {
        brdfs.add(reflection_oi);
        brdfs.add(transmission_oi);
      }

      /*! the ray transitions from inside to outside */
      else {
        brdfs.add(reflection_io);
        brdfs.add(transmission_io);
      }
    }

  private:
    DielectricReflection*   reflection_io;    //!< Reflection component for inside to outside transition.
    DielectricReflection*   reflection_oi;    //!< Reflection component for outside to inside transition.
    DielectricTransmission* transmission_io;  //!< Transmission component for inside to outside transition.
    DielectricTransmission* transmission_oi;  //!< Transmission component for outside to inside transition.
  };
}

#endif
