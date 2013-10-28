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
#include "brdfs/conductor.isph"
#include "brdfs/microfacet_metal.isph"

struct Metal
{
  Material base;
  float roughness;
  Conductor conductor;
  MicrofacetMetal microfacetMetal;
};

///////////////////////////////////////////////////////////////////////////////
// Implementation

void Metal__shade(const uniform Material *uniform _this,
                  const Ray&                  ray,       
                  const Medium&               currentMedium,
                  const DifferentialGeometry& dg,           
                  uniform CompositedBRDF&     brdfs)        
{
  const uniform Metal* uniform this = (const uniform Metal* uniform)_this;
  if (this->roughness == 0.0f)  CompositedBRDF__add(&brdfs,&this->conductor.base);
  else                          CompositedBRDF__add(&brdfs,&this->microfacetMetal.base);
}

void Metal__Constructor(uniform Metal *uniform this,
                        const uniform vec3f& reflectance,
                        const uniform vec3f& eta,
                        const uniform vec3f& k,
                        const uniform float roughness)
{
  Material__Constructor(&this->base,Material__Destructor,Metal__shade,NULL);
  this->roughness = roughness;
  Conductor__Constructor(&this->conductor,reflectance, eta, k);
  MicrofacetMetal__Constructor(&this->microfacetMetal,
                               reflectance, 
                               make_FresnelConductor(eta,k), 
                               make_PowerCosineDistribution(rcp(roughness)));
}

///////////////////////////////////////////////////////////////////////////////
// External API

export void* uniform Metal__new(const uniform vec3f& reflectance,
                                const uniform vec3f& eta,
                                const uniform vec3f& k,
                                const uniform float roughness)
{
  uniform Metal *uniform this = uniform new uniform Metal;
  Metal__Constructor(this,reflectance,eta,k,roughness);
  return this;
}
