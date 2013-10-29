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

#include "material.isph"
#include "brdfs/lambertian.isph"
#include "textures/texture.isph"

struct MatteTextured
{
  Material base;
  vec2f s0;                     //!< Offset for texture coordinates.
  vec2f ds;                     //!< Scaling for texture coordinates.
  uniform Texture* uniform Kd;  //!< Texture mapped to the surface.
};

///////////////////////////////////////////////////////////////////////////////
// Implementation

void MatteTextured__shade(const uniform Material *uniform _this,
                          const Ray&                  ray,    
                          const Medium&               currentMedium, 
                          const DifferentialGeometry& dg,            
                          uniform CompositedBRDF&     brdfs)         
{
  const uniform MatteTextured* uniform this = (const uniform MatteTextured* uniform) _this;
  COMPOSITED_BRDF_ADD(brdfs,Lambertian,
                      this->Kd->get(this->Kd,add(mul(this->ds,dg.st),this->s0)));
}

void MatteTextured__Destructor(uniform RefCount* uniform _this)
{
  uniform MatteTextured* uniform this = (uniform MatteTextured* uniform) _this;
  RefCount__DecRef(&this->Kd->base);
  Material__Destructor(_this);
}

inline void MatteTextured__Constructor(uniform MatteTextured *uniform this,
                                       uniform Texture* uniform Kd,
                                       const uniform vec2f& s0,
                                       const uniform vec2f& ds)
{
  Material__Constructor(&this->base,MatteTextured__Destructor,MatteTextured__shade,NULL);
  RefCount__IncRef(&Kd->base);
  this->Kd = Kd;
  this->s0 = s0;
  this->ds = ds;
}

///////////////////////////////////////////////////////////////////////////////
// External API

export void* uniform MatteTextured__new(void* uniform Kd,
                                        const uniform vec2f& s0,
                                        const uniform vec2f& ds)
{
  uniform MatteTextured *uniform this = uniform new uniform MatteTextured;
  MatteTextured__Constructor(this,(uniform Texture* uniform)Kd,s0,ds);
  return this;
}
