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

struct AmbientLight
{
  EnvironmentLight base;
  vec3f L;
};

AmbientLight* uniform AmbientLight__New(const uniform vec3f& L);

uniform Light* uniform AmbientLight__transform(const uniform Light *uniform _this, 
                                               const uniform AffineSpace3f& xfm) 
{
  const uniform AmbientLight* uniform this = (const uniform AmbientLight* uniform)_this;
  return (uniform Light* uniform) AmbientLight__New(this->L);
}

varying vec3f AmbientLight__Le(const uniform EnvironmentLight* uniform _this, varying const vec3f& wo) 
{
  const uniform AmbientLight* uniform this = (const uniform AmbientLight* uniform)_this;
  return this->L;
}

varying vec3f AmbientLight__eval(const uniform Light* uniform _this, 
                                 varying const DifferentialGeometry& dg, 
                                 varying const vec3f& wi) 
{
  const uniform AmbientLight* uniform this = (const uniform AmbientLight* uniform)_this;
  return this->L;
}

varying vec3f AmbientLight__sample(const uniform Light* uniform _this,
                                   varying const DifferentialGeometry& dg, 
                                   varying Sample3f& wi,
                                   varying float& tMax,
                                   varying const vec2f& s) 
{
  const uniform AmbientLight* uniform this = (const uniform AmbientLight* uniform)_this;
  wi = cosineSampleHemisphere(s.x, s.y, dg.Ns);
  tMax = 1e20f;
  return this->L;
}

void AmbientLight__Constructor(uniform AmbientLight* uniform this, const uniform vec3f L) 
{
  EnvironmentLight__Constructor(&this->base,Light__Destructor,ENV_LIGHT,
                                AmbientLight__transform,NULL,AmbientLight__eval,AmbientLight__sample,
                                AmbientLight__Le);
  this->L = L;
}

AmbientLight* uniform AmbientLight__New(const uniform vec3f& L)
{
  uniform AmbientLight* uniform light = uniform new uniform AmbientLight;
  AmbientLight__Constructor(light,L);
  return light;
}

export void* uniform AmbientLight__new(const uniform vec3f& L) {
  return AmbientLight__New(L);
}
