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
#include "brdfs/dielectriclayer_lambertian.isph"
#include "brdfs/dielectric_reflection.isph"
#include "brdfs/microfacet_plastic.isph"

struct Plastic
{
  Material base;
  float roughness;
  DielectricLayerLambertian substrate;
  DielectricReflection reflection;
  MicrofacetPlastic microfacet;
};

///////////////////////////////////////////////////////////////////////////////
// Implementation

void Plastic__shade(const uniform Material *uniform _this,
                  const Ray&                  ray,             
                  const Medium&               currentMedium,   
                  const DifferentialGeometry& dg,              
                  uniform CompositedBRDF&     brdfs)           
{
  const uniform Plastic* uniform this = (const uniform Plastic* uniform)_this;
  CompositedBRDF__add(&brdfs,&this->substrate.base);
  if (this->roughness == 0.0f)  CompositedBRDF__add(&brdfs,&this->reflection.base);
  else                          CompositedBRDF__add(&brdfs,&this->microfacet.base);
}

void Plastic__Constructor(uniform Plastic *uniform this,
                          const uniform vec3f& pigmentColor,
                          const uniform float eta,
                          const uniform float roughness)
{
  Material__Constructor(&this->base,Material__Destructor,Plastic__shade,NULL);
  this->roughness = roughness;

  DielectricLayerLambertian__Constructor(&this->substrate,
                                         make_vec3f(1.0f), 1.0f, eta, 
                                         make_Lambertian(pigmentColor));

  DielectricReflection__Constructor(&this->reflection, 1.0f, eta);

  MicrofacetPlastic__Constructor(&this->microfacet,make_vec3f(1.0f),
                                 make_FresnelDielectric(1.0f,eta), 
                                 make_PowerCosineDistribution(rcp(roughness)));
}

///////////////////////////////////////////////////////////////////////////////
// External API

export void* uniform Plastic__new(const uniform vec3f& pigmentColor,
                                  const uniform float eta,
                                  const uniform float roughness)
{
  uniform Plastic *uniform this = uniform new uniform Plastic;
  Plastic__Constructor(this,pigmentColor,eta,roughness);
  return this;
}
