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
#include "brdfs/dielectriclayer_lambertian.isph"
#include "brdfs/dielectriclayer_microfacet_metal.isph"

struct MetallicPaint
{
  Material base;
  DielectricReflection              reflection;
  DielectricLayerLambertian         paint;
  DielectricLayerMicrofacetMetal    glitter;
  bool                              hasGlitter;
};

///////////////////////////////////////////////////////////////////////////////
// Implementation

void MetallicPaint__shade(const uniform Material *uniform _this,
                          const Ray&                  ray,             
                          const Medium&               currentMedium,   
                          const DifferentialGeometry& dg,              
                          uniform CompositedBRDF&     brdfs)           
{
  uniform const MetallicPaint *uniform this = (uniform const MetallicPaint *uniform)_this;
  CompositedBRDF__add(&brdfs,&this->reflection.base);
  CompositedBRDF__add(&brdfs,&this->paint.base);
  if (this->hasGlitter) CompositedBRDF__add(&brdfs,&this->glitter.base);
}

inline void MetallicPaint__Constructor(uniform MetallicPaint *uniform this,
                                       const uniform vec3f& shadeColor,
                                       const uniform vec3f& glitterColor,
                                       const uniform float glitterSpread,
                                       const uniform float eta)
{
  Material__Constructor(&this->base,Material__Destructor,MetallicPaint__shade,NULL);
  DielectricReflection__Constructor(&this->reflection, 1.0f, eta);
  DielectricLayerLambertian__Constructor(&this->paint, make_vec3f(1.0f), 1.0f, eta, make_Lambertian(shadeColor));

  uniform vec3f etaAluminium = make_vec3f(0.62f,0.62f,0.62f);
  uniform vec3f kAluminium   = make_vec3f(4.8f,4.8f,4.8f);
  DielectricLayerMicrofacetMetal__Constructor(&this->glitter,
                                              make_vec3f(1.f), 1.f, eta, 
                                              make_MicrofacetMetal(glitterColor,
                                                                   make_FresnelConductor(etaAluminium,kAluminium),
                                                                   make_PowerCosineDistribution(rcp(glitterSpread))));
  this->hasGlitter = (glitterSpread != 0.0f & ne(glitterColor,make_vec3f(0.0f)));
}

///////////////////////////////////////////////////////////////////////////////
// External API

export void* uniform MetallicPaint__new(const uniform vec3f& shadeColor,
                                        const uniform vec3f& glitterColor,
                                        const uniform float glitterSpread,
                                        const uniform float eta)
{
  uniform MetallicPaint *uniform this = uniform new uniform MetallicPaint;
  MetallicPaint__Constructor(this,shadeColor,glitterColor,glitterSpread,eta);
  return this;
}
