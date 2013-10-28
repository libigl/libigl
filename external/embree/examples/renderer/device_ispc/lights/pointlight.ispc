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

struct PointLight
{
  uniform Light base;
  uniform vec3f P;
  uniform vec3f I;
};

PointLight* uniform PointLight__New(const uniform vec3f& P, const uniform vec3f& I);

uniform Light* uniform PointLight__transform(const uniform Light *uniform _this, 
                                             const uniform AffineSpace3f& xfm) 
{
  const uniform PointLight *uniform this = (const uniform PointLight *uniform)_this;
  return (uniform Light* uniform) PointLight__New(xfmPoint(xfm,this->P),this->I);
}

varying vec3f PointLight__eval(const uniform Light *uniform _this, 
                               varying const DifferentialGeometry &dg, 
                               varying const vec3f &wi) 
{
  return make_vec3f(0.0f);
}

varying vec3f PointLight__sample(const uniform Light *uniform _this,
                                   varying const DifferentialGeometry &dg, 
                                   varying Sample3f &wi,
                                   varying float &tMax,
                                   varying const vec2f &s) 
{
  const uniform PointLight *uniform this = (const uniform PointLight *uniform)_this;
  vec3f d = sub(this->P, dg.P);
  float distance = length(d);
  wi = make_Sample3f(mul(d,rcp(distance)), distance*distance);
  tMax = distance;
  return this->I;
}

void PointLight__Constructor(uniform PointLight* uniform this, const uniform vec3f& P, const uniform vec3f& I)
{
  Light__Constructor(&this->base,Light__Destructor,NORMAL_LIGHT,
                     PointLight__transform,NULL,PointLight__eval,PointLight__sample);
  this->P = P; 
  this->I = I;
}

PointLight* uniform PointLight__New(const uniform vec3f& P, const uniform vec3f& I)
{
  uniform PointLight *uniform this = uniform new uniform PointLight;
  PointLight__Constructor(this,P,I);
  return this;
}

export void* uniform PointLight__new(const uniform vec3f& P, const uniform vec3f& I) {
  return PointLight__New(P, I);
}
