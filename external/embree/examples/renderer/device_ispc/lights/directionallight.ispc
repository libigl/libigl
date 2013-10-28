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

#include "light.isph"
#include "samplers/shapesampler.isph"

struct DirectionalLight
{
  Light base;
  vec3f D;
  vec3f E;
};

DirectionalLight* uniform DirectionalLight__New(const uniform vec3f& D, const uniform vec3f& E);

uniform Light* uniform DirectionalLight__transform(const uniform Light *uniform _this, 
                                                   const uniform AffineSpace3f& xfm) 
{
  const uniform DirectionalLight *uniform this = (const uniform DirectionalLight *uniform)_this;
  return (uniform Light* uniform) DirectionalLight__New(xfmVector(xfm,this->D),this->E);
}

varying vec3f DirectionalLight__eval(const uniform Light *uniform _this, 
                                     varying const DifferentialGeometry &dg, 
                                     varying const vec3f &wi) {
  return make_vec3f(0.0f);
}

varying vec3f DirectionalLight__sample(const uniform Light *uniform _this,
                                       varying const DifferentialGeometry &dg, 
                                       varying Sample3f &wi,
                                       varying float &tMax,
                                       varying const vec2f &s) 
{
  const uniform DirectionalLight *uniform this = (const uniform DirectionalLight *uniform)_this;
  wi = make_Sample3f(neg(this->D),1.0f); 
  tMax = inf; 
  return this->E;
}

void DirectionalLight__Constructor (uniform DirectionalLight *uniform this,
                                    const uniform vec3f D,
                                    const uniform vec3f E)
{
  Light__Constructor(&this->base,Light__Destructor,NORMAL_LIGHT,
                     DirectionalLight__transform,NULL,DirectionalLight__eval,DirectionalLight__sample);
  this->D = normalize(D);
  this->E = E;
}

DirectionalLight* uniform DirectionalLight__New(const uniform vec3f& D, const uniform vec3f& E)
{
  uniform DirectionalLight *uniform this = uniform new uniform DirectionalLight;
  DirectionalLight__Constructor(this,D,E);
  return this;
}

export void* uniform DirectionalLight__new(const uniform vec3f& D, const uniform vec3f& E) {
  return DirectionalLight__New(D,E);
}
