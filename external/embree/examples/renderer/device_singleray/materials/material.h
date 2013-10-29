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

#ifndef __EMBREE_MATERIAL_H__
#define __EMBREE_MATERIAL_H__

#include "embree/common/ray.h"
#include "../shapes/differentialgeometry.h"
#include "../materials/medium.h"
#include "../brdfs/compositedbrdf.h"

namespace embree
{
  /*! Interface to different materials. A material implements a shader
   *  that creates a composition of BRDFs for the point to shade. */
  class Material : public RefCount {
    ALIGNED_CLASS
  public:

    /*! Material constructor. */
    Material (const Medium& mediumOutside = Medium::Vacuum(), /*!< The medium the geometry normal points into. */
              const Medium& mediumInside = Medium::Vacuum(),  /*!< The medium on the opposite side. */
              bool isMediaInterface = false)                  /*!< True if the surface is an interface between two different media. */
      : mediumOutside(mediumOutside), mediumInside(mediumInside), isMediaInterface(isMediaInterface) {}

    /*! Virtual destructor for materials */
    virtual ~Material() {}

    /*! Shades a location and and returns a composited BRDF */
    virtual void shade(const Ray&                  ray,              /*!< The ray arriving at the point to shade. */
                       const Medium&               currentMedium,    /*!< The medium this ray travels inside. */
                       const DifferentialGeometry& dg,               /*!< The point to shade on a surface. */
                       CompositedBRDF&             brdfs)            /*!< Container for generated BRDF components. */ const = 0;

    /*! Tracks the medium when crossing the surface. */
    __forceinline Medium nextMedium(const Medium& current) {
      if (!isMediaInterface) return current;
      return current == mediumInside ? mediumOutside : mediumInside;
    }

  public:
    Medium mediumOutside;   //!< Outside medium.
    Medium mediumInside;    //!< Inside medium.
    bool isMediaInterface;  //!< True if the surface is an interface between two different media.
  };
}

#endif
