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
#include "samplers/shapesampler.isph"

struct DepthOfFieldCamera
{
  Camera base;
  AffineSpace3f pixel2world;
  AffineSpace3f local2world;
  float lensRadius;
  float focalDistance;
};

///////////////////////////////////////////////////////////////////////////////
// Implementation

void DepthOfFieldCamera__ray(const uniform Camera *uniform _this,
                             varying Ray &ray,
                             const varying vec2f pixel, 
                             const varying vec2f sample)
{
  const uniform DepthOfFieldCamera *uniform this = (const uniform DepthOfFieldCamera *uniform)_this;
  const vec2f lens = uniformSampleDisk(sample,this->lensRadius);
  const vec3f begin = xfmPoint(this->local2world, make_vec3f(lens.x,lens.y,0.0f));
  const vec3f dir = add(mul(     pixel.x,this->pixel2world.l.vx), 
                        mul(1.0f-pixel.y,this->pixel2world.l.vy), 
                        this->pixel2world.l.vz);
  const vec3f end = add(this->pixel2world.p, mul(this->focalDistance,dir));
  init_Ray(ray,begin,normalize(sub(end,begin)));
}

void DepthOfFieldCamera__Constructor(uniform DepthOfFieldCamera* uniform this,
                                     const uniform AffineSpace3f local2world, 
                                     const uniform float& fov,
                                     const uniform float& aspectRatio,
                                     const uniform float& lensRadius,
                                     const uniform float& focalDistance)
{
  Camera__Constructor(&this->base,Camera__Destructor,DepthOfFieldCamera__ray);
  uniform vec3f W = xfmVector(local2world, make_vec3f(-0.5f*aspectRatio,-0.5f,0.5f*rcp(tan(deg2rad(0.5f*fov)))));
  uniform AffineSpace3f pixel2world = make_AffineSpace3f(mul(aspectRatio,local2world.l.vx),local2world.l.vy,W,local2world.p);
  this->pixel2world = pixel2world;
  this->local2world = local2world;
  this->lensRadius = lensRadius;
  this->focalDistance = focalDistance/length(add(mul(0.5f,pixel2world.l.vx), 
                                                 mul(0.5f,pixel2world.l.vy), 
                                                 pixel2world.l.vz));
}

///////////////////////////////////////////////////////////////////////////////
// External API

export void* uniform DepthOfFieldCamera__new(const uniform vec3f& vx, 
                                             const uniform vec3f& vy, 
                                             const uniform vec3f& vz, 
                                             const uniform vec3f& p, 
                                             const uniform float& fov,
                                             const uniform float& aspectRatio,
                                             const uniform float& lensRadius,
                                             const uniform float& focalDistance)
{
  uniform DepthOfFieldCamera *uniform this = uniform new uniform DepthOfFieldCamera;
  const uniform AffineSpace3f local2world = make_AffineSpace3f(vx,vy,vz,p);
  DepthOfFieldCamera__Constructor(this,local2world,fov,aspectRatio,lensRadius,focalDistance);
  return this;
}
