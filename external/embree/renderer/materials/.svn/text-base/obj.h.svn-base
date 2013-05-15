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

#ifndef __EMBREE_OBJ_H__
#define __EMBREE_OBJ_H__

#include "materials/material.h"
#include "brdfs/lambertian.h"
#include "brdfs/specular.h"
#include "brdfs/transmission.h"
#include "textures/texture.h"

namespace embree {

  /*! Material which models a subset of the features of the OBJ material format. */
  class Obj : public Material 
  {
  public:
    
    Obj (const Parms &options) 
    {
      /*! opacity texture and coefficient */
      map_d = options.getTexture("map_d");  d = options.getFloat("d", 1.0f);
      
      /*! diffuse reflectance texture and coefficient */
      map_Kd = options.getTexture("map_Kd");  Kd = options.getCol3f("Kd", one);
      
      /*! specular reflectance texture and coefficient */
      map_Ks = options.getTexture("map_Ks");  Ks = options.getCol3f("Ks", zero);
      
      /*! specular exponent texture and coefficient */
      map_Ns = options.getTexture("map_Ns");  Ns = options.getFloat("Ns", 10.0f);
    }
    
    void shade(const Ray &ray, const Medium &currentMedium, const DifferentialGeometry &dg, CompositedBRDF &brdfs) const 
    {
      /*! transmission */
      float d = this->d;  if (map_d) d *= map_d->get(dg.st).r; if (d < 1.0f) brdfs.add(NEW_BRDF(Transmission)(Col3f(1.0f - d)));
      
      /*! diffuse component */
      Col3f Kd = d*this->Kd;  if (map_Kd) Kd *= map_Kd->get(dg.st);  if (Kd != Col3f(zero)) brdfs.add(NEW_BRDF(Lambertian)(Kd));
      
      /*! specular exponent */
      float Ns = this->Ns;  if (map_Ns) Ns *= map_Ns->get(dg.st).r;
      
      /*! specular component */
      Col3f Ks = d*this->Ks;  if (map_Ks) Ks *= map_Ks->get(dg.st);  if (Ks != Col3f(zero)) brdfs.add(NEW_BRDF(Specular)(Ks, Ns));
    }
    
  private:
    Ref<Texture> map_d;   float d;     /*! opacity: 0 (transparent), 1 (opaque)                */
    Ref<Texture> map_Kd;  Col3f Kd;    /*! diffuse  reflectance: 0 (none), 1 (full)            */
    Ref<Texture> map_Ks;  Col3f Ks;    /*! specular reflectance: 0 (none), 1 (full)            */
    Ref<Texture> map_Ns;  float Ns;    /*! specular exponent: 0 (diffuse), infinity (specular) */
    Col3f Tf;                          /*! transmission: 0 (full absorption), 1 (transmissive) */
  };

}

#endif
