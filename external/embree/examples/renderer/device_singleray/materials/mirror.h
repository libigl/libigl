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

#ifndef __EMBREE_MIRROR_H__
#define __EMBREE_MIRROR_H__

#include "../materials/material.h"
#include "../brdfs/reflection.h"

namespace embree
{
  /*! Implements a mirror material. The reflected light can be
   *  modulated with a mirror reflectivity. */
  class Mirror : public Material
  {
  public:

    /*! Construction from parameters. */
    Mirror(const Parms& parms) {
      reflectance = parms.getColor("reflectance",one);
    }

    void shade(const Ray& ray, const Medium& currentMedium, const DifferentialGeometry& dg, CompositedBRDF& brdfs) const {
      brdfs.add(NEW_BRDF(Reflection)(reflectance));
    }

  protected:
    Color reflectance;  //!< Reflectivity of the mirror
  };
}

#endif
