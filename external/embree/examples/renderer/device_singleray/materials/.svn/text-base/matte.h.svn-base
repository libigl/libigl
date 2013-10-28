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

#ifndef __EMBREE_MATTE_H__
#define __EMBREE_MATTE_H__

#include "../materials/material.h"
#include "../brdfs/lambertian.h"

namespace embree
{
  /*! Implements a diffuse material. */
  class Matte : public Material
  {
  public:

    /*! Construction from parameters. */
    Matte (const Parms& parms) {
      reflectance = parms.getColor("reflectance",one);
    }

    void shade(const Ray& ray, const Medium& currentMedium, const DifferentialGeometry& dg, CompositedBRDF& brdfs) const {
      brdfs.add(NEW_BRDF(Lambertian)(reflectance));
    }

  protected:

    /*! Diffuse reflectance of the surface. The range is from 0
     *  (black) to 1 (white). */
    Color reflectance;
  };
}

#endif
