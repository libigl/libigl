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

struct Dielectric
{
  Material base;
  Medium mediumOutside;
  Medium mediumInside;
  DielectricReflection reflection_io;
  DielectricReflection reflection_oi;
  DielectricTransmission transmission_io;
  DielectricTransmission transmission_oi;
};

///////////////////////////////////////////////////////////////////////////////
// Implementation

void Dielectric__shade(const uniform Material *uniform _this,
                       const Ray&                  ray,      
                       const Medium&               currentMedium,
                       const DifferentialGeometry& dg,           
                       uniform CompositedBRDF&     brdfs)        
{
  uniform const Dielectric *uniform this = (uniform const Dielectric *uniform)_this;

  /*! the ray transitions from outside to inside */
  if (eq(currentMedium,this->mediumOutside)) {
    CompositedBRDF__add(&brdfs,&this->reflection_oi.base);
    CompositedBRDF__add(&brdfs,&this->transmission_oi.base);
  }
  
  /*! the ray transitions from inside to outside */
  else {
    CompositedBRDF__add(&brdfs,&this->reflection_io.base);
    CompositedBRDF__add(&brdfs,&this->transmission_io.base);
  }
}

void Dielectric__selectNextMedium(const uniform Material *uniform _this,
                                  varying Medium &lastMedium) 
{
  const uniform Dielectric *uniform this = (const uniform Dielectric *uniform)_this;
  
  if (eq(lastMedium,this->mediumOutside)) lastMedium = this->mediumInside;
  else                                    lastMedium = this->mediumOutside;
}

void Dielectric__Constructor(uniform Dielectric* uniform this,
                             const uniform float etaOutside,
                             const uniform vec3f& transmissionOutside,
                             const uniform float etaInside,
                             const uniform vec3f& transmissionInside) 
{
  Material__Constructor(&this->base,Material__Destructor,Dielectric__shade,Dielectric__selectNextMedium);
  this->mediumOutside.eta          = etaOutside;
  this->mediumOutside.transmission = transmissionOutside;
  this->mediumInside.eta           = etaInside;
  this->mediumInside.transmission  = transmissionInside;
  DielectricReflection__Constructor  (&this->reflection_io  ,etaInside ,etaOutside);
  DielectricReflection__Constructor  (&this->reflection_oi  ,etaOutside,etaInside);
  DielectricTransmission__Constructor(&this->transmission_io,etaInside ,etaOutside);
  DielectricTransmission__Constructor(&this->transmission_oi,etaOutside,etaInside);
}

///////////////////////////////////////////////////////////////////////////////
// External API

export void* uniform Dielectric__new(const uniform float etaOutside,
                                     const uniform vec3f& transmissionOutside,
                                     const uniform float etaInside,
                                     const uniform vec3f& transmissionInside) 
{
  uniform Dielectric *uniform mat = uniform new uniform Dielectric;
  Dielectric__Constructor(mat,etaOutside,transmissionOutside,etaInside,transmissionInside);
  return mat;
}
