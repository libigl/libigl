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

#include "tonemapper.isph"

struct DefaultToneMapper
{
  ToneMapper base;
  float gamma;
  float rcpGamma;
  bool vignetting;
};

varying vec3f DefaultToneMapper__ToneMap(const uniform ToneMapper *uniform _this, 
                                         const varying vec3f v, 
                                         const varying int x, const varying int y, const uniform vec2ui size) 
{ 
  const uniform DefaultToneMapper *uniform this = (const uniform DefaultToneMapper *uniform)_this;
  
  vec3f c = v;
  if (this->gamma != 1.0f) 
    c = pow(c,this->rcpGamma);
  if (this->vignetting) {
    const vec2f center = make_vec2f(0.5f*(float)size.x,0.5f*(float)size.y);
    float d = length(sub(make_vec2f((float)x,(float)y),center)) / (size.x*0.5f);
    //c *= pow(cos(d*0.5f),3.0f);
    c = mul(c,pow(cos(d*0.5f),3.0f));
  }
  return c;
}

void DefaultToneMapper__Constructor(uniform DefaultToneMapper* uniform this, uniform float gamma, uniform bool vignetting)
{
  ToneMapper__Constructor(&this->base,ToneMapper__Destructor,DefaultToneMapper__ToneMap);
  this->rcpGamma = rcp(gamma);
  this->vignetting = vignetting;
}

export void* uniform DefaultToneMapper__new(const uniform float& gamma, const uniform bool& vignetting)
{
  uniform DefaultToneMapper *uniform obj = uniform new uniform DefaultToneMapper;
  DefaultToneMapper__Constructor(obj,gamma,vignetting);
  return obj;
}
