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
#include "brdfs/minneart.isph"
#include "brdfs/velvety.isph"

struct Velvet 
{
  uniform Material base;

  /*! Diffuse reflectance of the surface. The range is from 0
   *  (black) to 1 (white). */
  vec3f reflectance;
  
  /*! Amount of back scattering. The range is from 0 (no back
   *  scattering) to inf (maximum back scattering). */
  float backScattering;
  
  /*! Color of horizon scattering. */
  vec3f horizonScatteringColor;
  
  /*! Fall-off of horizon scattering. */
  float horizonScatteringFallOff;
  
  // the two BRDFs
  Minneart minneart;
  Velvety  velvety;
};

///////////////////////////////////////////////////////////////////////////////
// Implementation

void Velvet__shade(const uniform Material *uniform _this,
                  const Ray&                  ray,            
                  const Medium&               currentMedium,  
                  const DifferentialGeometry& dg,             
                  uniform CompositedBRDF&     brdfs)          
{
  const uniform Velvet* uniform this = (const uniform Velvet* uniform)_this;
  CompositedBRDF__add(&brdfs,&this->minneart.base);
  CompositedBRDF__add(&brdfs,&this->velvety.base);
}

void Velvet__Constructor(uniform Velvet *uniform this,
                        const uniform vec3f& reflectance,
                        const uniform vec3f& horizonScatteringColor,
                        const uniform float horizonScatteringFallOff,
                        const uniform float backScattering)
{
  Material__Constructor(&this->base,Material__Destructor,Velvet__shade,NULL);
  this->reflectance = reflectance;
  this->horizonScatteringColor = horizonScatteringColor;
  this->horizonScatteringFallOff = horizonScatteringFallOff;
  this->backScattering = backScattering;
  Minneart__Constructor(&this->minneart,reflectance,backScattering);
  Velvety__Constructor (&this->velvety,horizonScatteringColor,horizonScatteringFallOff);
}

///////////////////////////////////////////////////////////////////////////////
// External API

export void* uniform Velvet__new(const uniform vec3f& reflectance,
                                 const uniform vec3f& horizonScatteringColor,
                                 const uniform float horizonScatteringFallOff,
                                 const uniform float backScattering)
{
  uniform Velvet *uniform this = uniform new uniform Velvet;
  Velvet__Constructor(this,reflectance,horizonScatteringColor,horizonScatteringFallOff,backScattering);
  return this;
}
