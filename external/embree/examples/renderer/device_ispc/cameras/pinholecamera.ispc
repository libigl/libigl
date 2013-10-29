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

#include "camera.isph"

struct PinHoleCamera
{
  Camera base;
  AffineSpace3f pixel2world;
};

///////////////////////////////////////////////////////////////////////////////
// Implementation

void PinHoleCamera__ray(const uniform Camera *uniform _this,	
                        varying Ray &ray,	      
                        const varying vec2f pixel,
                        const varying vec2f sample)
{
  const uniform PinHoleCamera *uniform this = (const uniform PinHoleCamera *uniform)_this;
  init_Ray(ray,this->pixel2world.p,normalize(add(mul(pixel.x,this->pixel2world.l.vx), 
                                                 mul(1.0f-pixel.y,this->pixel2world.l.vy), 
                                                 this->pixel2world.l.vz)));
}

void PinHoleCamera__Constructor(uniform PinHoleCamera *uniform this,
                                const uniform AffineSpace3f local2world, 
                                const uniform float& fov,
                                const uniform float& aspectRatio)
{
  Camera__Constructor(&this->base,Camera__Destructor,PinHoleCamera__ray);
  uniform vec3f W = xfmVector(local2world, make_vec3f(-0.5f*aspectRatio,-0.5f,0.5f*rcp(tan(deg2rad(0.5f*fov)))));
  this->pixel2world = make_AffineSpace3f(mul(aspectRatio,local2world.l.vx),local2world.l.vy,W,local2world.p);
}

///////////////////////////////////////////////////////////////////////////////
// External API

export void* uniform PinHoleCamera__new(const uniform vec3f& vx, 
                                        const uniform vec3f& vy, 
                                        const uniform vec3f& vz, 
                                        const uniform vec3f& p, 
                                        const uniform float& fov,
                                        const uniform float& aspectRatio)
{
  uniform PinHoleCamera *uniform this = uniform new uniform PinHoleCamera;
  const uniform AffineSpace3f local2world = make_AffineSpace3f(vx,vy,vz,p);
  PinHoleCamera__Constructor(this,local2world,fov,aspectRatio);
  return this;
}
