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

#ifndef __EMBREE_METALLIC_PAINT_H__
#define __EMBREE_METALLIC_PAINT_H__

#include "../materials/material.h"
#include "../brdfs/microfacet.h"
#include "../brdfs/dielectric.h"
#include "../brdfs/dielectriclayer.h"
#include "../brdfs/lambertian.h"

namespace embree
{
  /*! Implements a car paint BRDF. Models a dielectric layer over a
   *  diffuse ground layer. Additionally the ground layer may contain
   *  metallic glitter. */
  class MetallicPaint : public Material
  {
    typedef Microfacet<FresnelConductor,PowerCosineDistribution> MicrofacetGlitter;

  public:

    /*! Construction from parameters. */
    MetallicPaint (const Parms& parms)
    {
      /*! extract parameters */
      shadeColor    = parms.getColor("shadeColor",one);
      glitterColor  = parms.getColor("glitterColor",zero);
      glitterSpread = parms.getFloat("glitterSpread",1.0f);
      eta           = parms.getFloat("eta",1.4f);

      /*! precompute BRDF component for performance reasons */
      reflection = new DielectricReflection(1.0f, eta);
      paint      = new DielectricLayer<Lambertian >(one, 1.0f, eta, Lambertian(shadeColor));
    }

    /*! Destructor */
    ~MetallicPaint ()
    {
      if (reflection) delete reflection; reflection = NULL;
      if (paint     ) delete paint     ; paint      = NULL;
    }
    
    void shade(const Ray& ray, const Medium& currentMedium, const DifferentialGeometry& dg, CompositedBRDF& brdfs) const
    {
      brdfs.add(reflection);
      brdfs.add(paint);

      if (glitterSpread != 0 && glitterColor != Color(zero)) 
      {
        /*! Use the fresnel relfectance of Aluminium for the flakes. Modulate with glitterColor. */
        Color etaAluminium(0.62f,0.62f,0.62f);
        Color kAluminium(4.8,4.8,4.8);
        brdfs.add(NEW_BRDF(DielectricLayer<MicrofacetGlitter>)(one, 1.0f, eta, MicrofacetGlitter(glitterColor,
                                                                                                  FresnelConductor(etaAluminium,kAluminium),
                                                                                                  PowerCosineDistribution(rcp(glitterSpread),dg.Ns))));
      }
    }

  protected:
    Color glitterColor;
    float glitterSpread;
    Color shadeColor;
    float eta;
    DielectricReflection*        reflection;  //!< Precomputed dielectric reflection component.
    DielectricLayer<Lambertian>* paint;       //!< Diffuse layer covered by dielectricum.
  };
}

#endif
