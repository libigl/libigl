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

#ifndef __EMBREE_OBJ_H__
#define __EMBREE_OBJ_H__

#include "../materials/material.h"
#include "../brdfs/lambertian.h"
#include "../brdfs/specular.h"
#include "../brdfs/transmission.h"
#include "../textures/texture.h"

namespace embree {

    /*! material which models a subset of the features of the OBJ material format */
    class Obj : public Material {
    public:

        Obj(const Parms &options) {

            /*! opacity texture and coefficient */
            map_d = options.getTexture("map_d");  d = options.getFloat("d", 1.0f);

            /*! diffuse reflectance texture and coefficient */
            map_Kd = options.getTexture("map_Kd");  Kd = options.getColor("Kd", one);

            /*! specular reflectance texture and coefficient */
            map_Ks = options.getTexture("map_Ks");  Ks = options.getColor("Ks", zero);

            /*! specular exponent texture and coefficient */
            map_Ns = options.getTexture("map_Ns");  Ns = options.getFloat("Ns", 10.0f);

            /*! get bump map */
            map_Bump = options.getTexture("map_Bump");
        }

        void shade(const Ray &ray, const Medium &currentMedium, const DifferentialGeometry &dg, CompositedBRDF &brdfs) const 
        {
          if (this->map_Bump) {
            const Color bump = map_Bump->get(dg.st);
            const Vector3f b(2.0f*bump.r-1.0f,2.0f*bump.g-1.0f,2.0f*bump.b-1.0f);
            dg.Ns = normalize(b.x*dg.Tx + b.y*dg.Ty + b.z*dg.Ns);
          }

          /*! transmission */
          float d = this->d;  if (map_d) d *= map_d->get(dg.st).r; if (d < 1.0f) brdfs.add(NEW_BRDF(Transmission)(Color(1.0f - d)));
          
          /*! diffuse component */
          Color Kd = d*this->Kd;  if (map_Kd) Kd *= map_Kd->get(dg.st);  if (Kd != Color(zero)) brdfs.add(NEW_BRDF(Lambertian)(Kd));
          
          /*! specular exponent */
          float Ns = this->Ns;  if (map_Ns) Ns *= map_Ns->get(dg.st).r;
          
          /*! specular component */
          Color Ks = d*this->Ks;  if (map_Ks) Ks *= map_Ks->get(dg.st);  if (Ks != Color(zero)) brdfs.add(NEW_BRDF(Specular)(Ks, Ns));
        }

    protected:
        Ref<Texture> map_d;   float d;     /*! opacity: 0 (transparent), 1 (opaque)                */
        Ref<Texture> map_Kd;  Color Kd;    /*! diffuse  reflectance: 0 (none), 1 (full)            */
        Ref<Texture> map_Ks;  Color Ks;    /*! specular reflectance: 0 (none), 1 (full)            */
        Ref<Texture> map_Ns;  float Ns;    /*! specular exponent: 0 (diffuse), infinity (specular) */
        Ref<Texture> map_Bump;             /*! bump map */
    };
}

#endif // __EMBREE_OBJ_H__

