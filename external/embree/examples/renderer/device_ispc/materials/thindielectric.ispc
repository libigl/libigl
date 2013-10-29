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
#include "brdfs/dielectric_reflection.isph"
#include "brdfs/dielectric_transmission.isph"

struct ThinDielectric
{
  uniform Material base;
  DielectricReflection reflection;
  ThinDielectricTransmission transmission;
};

///////////////////////////////////////////////////////////////////////////////
// Implementation

void ThinDielectric__shade(const uniform Material *uniform _this,
                       const Ray&                  ray,            
                       const Medium&               currentMedium,  
                       const DifferentialGeometry& dg,             
                       uniform CompositedBRDF&     brdfs)          
{
  uniform const ThinDielectric *uniform this = (uniform const ThinDielectric *uniform)_this;
  CompositedBRDF__add(&brdfs,&this->reflection.base);
  CompositedBRDF__add(&brdfs,&this->transmission.base);
}

void ThinDielectric__Constructor(uniform ThinDielectric* uniform this,
                                 const uniform vec3f& transmission,
                                 const uniform float eta,
                                 const uniform float thickness)
{
  Material__Constructor(&this->base,Material__Destructor,ThinDielectric__shade,NULL);
  DielectricReflection__Constructor      (&this->reflection  , 1.0f, eta);
  ThinDielectricTransmission__Constructor(&this->transmission, 1.0f, eta, transmission, thickness);
}

///////////////////////////////////////////////////////////////////////////////
// External API

export void* uniform ThinDielectric__new(const uniform vec3f& transmission,
                                         const uniform float eta,
                                         const uniform float thickness) 
{
  uniform ThinDielectric *uniform this = uniform new uniform ThinDielectric;
  ThinDielectric__Constructor(this,transmission,eta,thickness);
  return this;
}
