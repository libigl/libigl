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

struct SpotLight
{
  Light base;
  vec3f P;             //!< Position of the spot light
  vec3f D;             //!< Light direction of the spot light
  vec3f I;             //!< Radiant intensity (W/sr)
  float angleMin;      //!< Start of linear falloff region
  float angleMax;      //!< End of linear falloff region
  float cosAngleMin;   //!< Start of linear falloff region
  float cosAngleMax;   //!< End of linear falloff region
};

SpotLight* uniform SpotLight__New(const uniform vec3f& P,
                                  const uniform vec3f& D,
                                  const uniform vec3f& I,
                                  const uniform float angleMin,
                                  const uniform float angleMax);

uniform Light* uniform SpotLight__transform(const uniform Light *uniform _this, 
                                            const uniform AffineSpace3f& xfm) 
{
  const uniform SpotLight *uniform this = (const uniform SpotLight *uniform)_this;
  return (uniform Light* uniform) SpotLight__New(xfmPoint(xfm,this->P),xfmVector(xfm,this->D),this->I,this->angleMin,this->angleMax);
}

varying vec3f SpotLight__eval(const uniform Light *uniform _this, 
                              varying const DifferentialGeometry &dg, 
                              varying const vec3f &wi) 
{
  return make_vec3f(0.0f);
}

varying vec3f SpotLight__sample(const uniform Light *uniform _this,
                                varying const DifferentialGeometry &dg, 
                                varying Sample3f &wi,
                                varying float &tMax,
                                varying const vec2f &s) 
{
  const uniform SpotLight *uniform this = (const uniform SpotLight *uniform)_this;
  const vec3f d = sub(this->P, dg.P);
  const float distance = length(d);
  wi = make_Sample3f(mul(d,rcp(distance)), distance*distance);
  tMax = distance;

  const float cosAngle = -dot(wi.v,this->D);
  if (this->cosAngleMin != this->cosAngleMax)
    return mul(this->I, clamp((cosAngle - this->cosAngleMax)*rcp(this->cosAngleMin - this->cosAngleMax)));
  else if (cosAngle > this->cosAngleMin)
    return this->I;
  else 
    return make_vec3f(0.0f);
}

void SpotLight__Constructor (uniform SpotLight *uniform this,
                             const uniform vec3f P,
                             const uniform vec3f D,
                             const uniform vec3f I,
                             const uniform float angleMin,
                             const uniform float angleMax)
{
  Light__Constructor(&this->base,Light__Destructor,NORMAL_LIGHT,
                     SpotLight__transform,NULL,SpotLight__eval,SpotLight__sample);
  this->P = P;
  this->D = normalize(D);
  this->I = I;
  this->angleMin = angleMin;
  this->angleMax = angleMax;
  this->cosAngleMin = cos(0.5f*deg2rad(angleMin));
  this->cosAngleMax = cos(0.5f*deg2rad(angleMax));
}

SpotLight* uniform SpotLight__New(const uniform vec3f& P,
                                  const uniform vec3f& D,
                                  const uniform vec3f& I,
                                  const uniform float angleMin,
                                  const uniform float angleMax)
{
  uniform SpotLight *uniform this = uniform new uniform SpotLight;
  SpotLight__Constructor(this,P,D,I,angleMin,angleMax);
  return this;
}

export void* uniform SpotLight__new(const uniform vec3f& P,
                                    const uniform vec3f& D,
                                    const uniform vec3f& I,
                                    const uniform float angleMin,
                                    const uniform float angleMax)
{
  return SpotLight__New(P,D,I,angleMin,angleMax);
}
