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

struct DistantLight
{
  EnvironmentLight base;
  vec3f D;             //!< Light direction
  vec3f L;             //!< Radiance (W/(m^2*sr))
  float halfAngle;     //!< Half illumination angle
  float radHalfAngle;  //!< Half illumination angle in radians
  float cosHalfAngle;  //!< Cosine of half illumination angle
};

DistantLight* uniform DistantLight__New(const uniform vec3f& D,
                                        const uniform vec3f& L,
                                        const uniform float halfAngle);

uniform Light* uniform DistantLight__transform(const uniform Light *uniform _this, 
                                               const uniform AffineSpace3f& xfm) 
{
  const uniform DistantLight *uniform this = (const uniform DistantLight *uniform)_this;
  return (uniform Light* uniform) DistantLight__New(xfmVector(xfm,this->D),this->L,this->halfAngle);
}

varying vec3f DistantLight__Le(const uniform EnvironmentLight *uniform _this, varying const vec3f &wo) 
{
  const uniform DistantLight *uniform this = (const uniform DistantLight *uniform)_this;
  if (dot(wo,this->D) >= this->cosHalfAngle) return this->L;
  return make_vec3f(0.0f);
}

varying vec3f DistantLight__eval(const uniform Light *uniform _this, 
                                 varying const DifferentialGeometry &dg, 
                                 varying const vec3f &wi)
{
  const uniform DistantLight *uniform this = (const uniform DistantLight *uniform)_this;
  if (-dot(wi,this->D) >= this->cosHalfAngle) return this->L;
  return make_vec3f(0.0f);
}

varying vec3f DistantLight__sample(const uniform Light *uniform _this,
                                   varying const DifferentialGeometry &dg, 
                                   varying Sample3f &wi,
                                   varying float &tMax,
                                   varying const vec2f &s) 
{
  const uniform DistantLight *uniform this = (const uniform DistantLight *uniform)_this;
  wi = uniformSampleCone(s.x,s.y,this->radHalfAngle,neg(this->D));
  tMax = 1e20f;
  return this->L;
}

void DistantLight__Constructor(uniform DistantLight *uniform this,
                               const uniform vec3f D,
                               const uniform vec3f L,
                               const uniform float halfAngle)
{
  EnvironmentLight__Constructor(&this->base,Light__Destructor,ENV_LIGHT,
                                DistantLight__transform,NULL,DistantLight__eval,DistantLight__sample,
                                DistantLight__Le);
  this->D = normalize(D);
  this->L = L;
  this->halfAngle = halfAngle;
  this->radHalfAngle = deg2rad(halfAngle);
  this->cosHalfAngle = cos(deg2rad(halfAngle));
}

DistantLight* uniform DistantLight__New(const uniform vec3f& D,
                                        const uniform vec3f& L,
                                        const uniform float halfAngle)
{
  uniform DistantLight *uniform this = uniform new uniform DistantLight;
  DistantLight__Constructor(this,D,L,halfAngle);
  return this;
}

export void* uniform DistantLight__new(const uniform vec3f& D,
                                       const uniform vec3f& L,
                                       const uniform float halfAngle)
{
  return DistantLight__New(D,L,halfAngle);
}
