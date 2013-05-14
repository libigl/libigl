// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#include "materials/material.h"
#include "brdfs/minnaert.h"
#include "brdfs/velvety.h"

namespace embree
{
  /*! Implements a velvet material. */
  class Velvet : public Material
  {
  public:

    /*! Construction from parameters. */
    Velvet (const Parms& parms) {
      reflectance = parms.getCol3f("reflectance",one);
      backScattering = parms.getFloat("backScattering",zero);
      horizonScatteringColor = parms.getCol3f("horizonScatteringColor",one);
      horizonScatteringFallOff = parms.getFloat("horizonScatteringFallOff",zero);
    }

    void shade(const Ray& ray, const Medium& currentMedium, const DifferentialGeometry& dg, CompositedBRDF& brdfs) const {
      brdfs.add(NEW_BRDF(Minnaert)(reflectance, backScattering));
      brdfs.add(NEW_BRDF(Velvety)(horizonScatteringColor, horizonScatteringFallOff));
    }

  private:

    /*! Diffuse reflectance of the surface. The range is from 0
     *  (black) to 1 (white). */
    Col3f reflectance;

    /*! Amount of back scattering. The range is from 0 (no back
     *  scattering) to inf (maximum back scattering). */
    float backScattering;

    /*! Color of horizon scattering. */
    Col3f horizonScatteringColor;

    /*! Fall-off of horizon scattering. */
    float horizonScatteringFallOff;
  };

  /*! Implements a velvet material. */
  class Velvet2 : public Material
  {
  public:

    /*! Construction from parameters. */
    Velvet2 (const Parms& parms) {
      contrast = parms.getFloat("contrast",0.5f);
      s0 = parms.getVec2f("s0",Vec2f(0.0f,0.0f));
      ds = parms.getVec2f("ds",Vec2f(1.0f,1.0f));
      Kd = parms.getTexture("Kd");
      reflectance = parms.getCol3f("reflectance",one);
      backScattering = parms.getFloat("backScattering",zero);
      horizonScatteringColor = parms.getCol3f("horizonScatteringColor",one);
      horizonScatteringFallOff = parms.getFloat("horizonScatteringFallOff",zero);
    }

    void shade(const Ray& ray, const Medium& currentMedium, const DifferentialGeometry& dg, CompositedBRDF& brdfs) const {
      const Col3f color = contrast*(Kd->get(s0+ds*dg.st)-Col3f(0.5f))+Col3f(0.5f);
      brdfs.add(NEW_BRDF(Minnaert)(reflectance*color, backScattering));
      brdfs.add(NEW_BRDF(Velvety)(horizonScatteringColor*color, horizonScatteringFallOff));
    }

  private:

    /*! Diffuse reflectance of the surface. The range is from 0
     *  (black) to 1 (white). */
    Col3f reflectance;

    /*! Amount of back scattering. The range is from 0 (no back
     *  scattering) to inf (maximum back scattering). */
    float backScattering;

    /*! Color of horizon scattering. */
    Col3f horizonScatteringColor;

    /*! Fall-off of horizon scattering. */
    float horizonScatteringFallOff;

    float contrast;
    Vec2f s0;         //!< Offset for texture coordinates.
    Vec2f ds;         //!< Scaling for texture coordinates.
    Ref<Texture> Kd;
  };
}

#endif
