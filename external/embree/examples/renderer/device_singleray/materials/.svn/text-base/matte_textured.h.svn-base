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

#ifndef __EMBREE_MATTE_TEXTURED_H__
#define __EMBREE_MATTE_TEXTURED_H__

#include "../materials/material.h"
#include "../brdfs/lambertian.h"
#include "../textures/texture.h"

namespace embree
{
  /*! Implements a diffuse and textured material.*/
  class MatteTextured : public Material
  {
  public:

    /*! Construction from parameters. */
    MatteTextured (const Parms& parms)
    {
      Kd = parms.getTexture("Kd");
      s0 = parms.getVec2f("s0",Vec2f(0.0f,0.0f));
      ds = parms.getVec2f("ds",Vec2f(1.0f,1.0f));
    }

    void shade(const Ray& ray, const Medium& currentMedium, const DifferentialGeometry& dg, CompositedBRDF& brdfs) const {
      if (Kd) brdfs.add(NEW_BRDF(Lambertian)(Kd->get(ds*dg.st+s0)));
    }

  protected:
    Vec2f s0;         //!< Offset for texture coordinates.
    Vec2f ds;         //!< Scaling for texture coordinates.
    Ref<Texture> Kd;  //!< Texture mapped to the surface.
  };
}

#endif
