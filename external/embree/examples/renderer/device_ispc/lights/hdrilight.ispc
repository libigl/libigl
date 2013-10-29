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
#include "textures/image.isph"
#include "samplers/distribution2d.isph"

struct HDRILight
{
  EnvironmentLight base;
  AffineSpace3f local2world;          //!< Transformation from light space into world space
  AffineSpace3f world2local;          //!< Transformation from world space into light space
  unsigned int width, height;         //!< Width and height of the used image.
  vec3f L;                            //!< Scaling factor for the image.
  uniform Image *image;              //!< The image mapped to the environment.
  Distribution2D* uniform distribution;        //!< The 2D distribution used to importance sample.
};

void HDRILight__Constructor(uniform HDRILight *uniform this,
                            const uniform AffineSpace3f& local2world,
                            const uniform vec3f& L,
                            uniform Image *uniform image,
                            uniform Distribution2D* uniform distribution);

uniform Light* uniform HDRILight__transform(const uniform Light *uniform _this, 
                                            const uniform AffineSpace3f& xfm) 
{
  const uniform HDRILight *uniform this = (const uniform HDRILight *uniform)_this;
  uniform HDRILight *uniform light = uniform new uniform HDRILight;
  HDRILight__Constructor(light,mul(xfm,this->local2world),this->L,this->image,this->distribution);
  return &light->base.base;
}

varying vec3f HDRILight__Le(const uniform EnvironmentLight *uniform _this, const varying vec3f &wo) 
{
  const uniform HDRILight *uniform this = (const uniform HDRILight *uniform)_this;

  const vec3f wi = xfmVector(this->world2local, neg(wo));

  const float theta = acos(clamp(wi.y,-1.0f,1.0f));
  float phi = atan2(-wi.z,-wi.x);
  if (phi < 0.f) phi += 2.0f * (float)(M_PI);
  const float u = 1.0f - (phi * (float)(one_over_two_pi));
  const float v = theta * (float)(one_over_pi);

  return mul(this->L, this->image->get_bilinear_varying(this->image,u,v));
}

varying vec3f HDRILight__eval(const uniform Light *uniform _this, 
                               varying const DifferentialGeometry &dg, 
                               varying const vec3f &wi) 
{
  const uniform HDRILight *uniform this = (const uniform HDRILight *uniform)_this;
  return HDRILight__Le(&this->base,neg(wi));
}

varying vec3f HDRILight__sample(const uniform Light *uniform _this,
                                varying const DifferentialGeometry &dg, 
                                varying Sample3f &wi,
                                varying float &tMax,
                                varying const vec2f &sample)
{
  const uniform HDRILight *uniform this = (const uniform HDRILight *uniform)_this;

  const Sample2f pixelF = Distribution2D__sample(this->distribution,sample);
  const float theta = (float)(M_PI) * pixelF.v.y*rcp((float)(this->height));
  const float phi   = (float)(two_pi) * (1.0f - pixelF.v.x*rcp((float)(this->width)));
  const float sin_theta = sin(theta);
  const float neg_sin_theta = -sin_theta;
  const vec3f _wi = make_vec3f(neg_sin_theta*cos(phi),
                               cos(theta),neg_sin_theta*sin(phi));
  wi = make_Sample3f(xfmVector(this->local2world, _wi),
                     pixelF.pdf*rcp(((float)(two_pi) * (float)(pi)) * sin_theta));
  tMax = inf;

  return mul(this->L, this->image->get_nearest_varying(this->image,
                                                       clamp((int)(pixelF.v.x), 0, (int)(this->width-1)),
                                                       clamp((int)(pixelF.v.y), 0, (int)(this->height-1))));
}

void HDRILight__Destructor(uniform RefCount* uniform _this)
{ 
  uniform HDRILight* uniform this = (uniform HDRILight* uniform) _this;
  RefCount__DecRef(&this->image->base);
  RefCount__DecRef(&this->distribution->base);
  Light__Destructor(_this);
}

void HDRILight__Constructor(uniform HDRILight *uniform this,
                            const uniform AffineSpace3f& local2world,
                            const uniform vec3f& L,
                            uniform Image *uniform image,
                            uniform Distribution2D* uniform distribution)
{
  EnvironmentLight__Constructor(&this->base,HDRILight__Destructor,
                                (LightType)(ENV_LIGHT | PRECOMPUTED_LIGHT),
                                HDRILight__transform,NULL,HDRILight__eval,HDRILight__sample,HDRILight__Le);

  RefCount__IncRef(&image->base);
  this->image  = image;
  this->width  = image->size.x;
  this->height = image->size.y;

  RefCount__IncRef(&distribution->base);
  this->distribution = distribution;
  
  this->L      = L;
  this->local2world = local2world;
  this->world2local = rcp(local2world);
}

void HDRILight__Constructor(uniform HDRILight *uniform this,
                            const uniform AffineSpace3f &local2world,
                            const uniform vec3f &L,
                            uniform Image* uniform image)
{
  const uniform uint width  = image->size.x;
  const uniform uint height = image->size.y;

  /* calculate importance */
  uniform float* uniform importance = uniform new uniform float[width*height];  
  for (uniform int y = 0; y < height; y++) {
    uniform float tmp = sin(pi*(y+0.5f)*rcp((float)height));
    for (uniform int x = 0; x < width; x++) {
      importance[y*width+x] = tmp * reduce_add(image->get_nearest_uniform(image,x,y));
    }
  }

  /* create distribution */
  uniform Distribution2D* uniform distribution = Distribution2D__new(importance,image->size);
  delete[] importance;

  /* call constructor */
  HDRILight__Constructor(this,local2world,L,image,distribution);
}

export void* uniform HDRILight__new(const uniform vec3f& vx,
                                    const uniform vec3f& vy,
                                    const uniform vec3f& vz,
                                    const uniform vec3f& p,
                                    const uniform vec3f& L,
                                    void* uniform image)
{
  const uniform AffineSpace3f local2world = make_AffineSpace3f(vx,vy,vz,p);
  uniform HDRILight *uniform this = uniform new uniform HDRILight;
  HDRILight__Constructor(this,local2world,L,(uniform Image *uniform)image);
  return this;
}
