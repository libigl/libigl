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

#include "hdrilight.h"
#include "sys/stl/array2d.h"

namespace embree
{
  HDRILight::HDRILight(const AffineSpace3f& local2world, 
                       unsigned width, unsigned height, 
                       const Color& L, const Ref<Image>& pixels, 
                       const Ref<Distribution2D>& distribution,
                       light_mask_t illumMask,
                       light_mask_t shadowMask)
    : EnvironmentLight(illumMask,shadowMask), 
      local2world(local2world), 
      world2local(rcp(local2world)), 
      width(width), 
      height(height), 
      L(L), 
      pixels(pixels), 
      distribution(distribution) 
  {}
  
  HDRILight::HDRILight(const Parms& parms)
    : width(0), height(0), pixels(null)
  {
    local2world = parms.getTransform("local2world",one);
    world2local = rcp(local2world);
    L = parms.getColor("L",one);
    pixels = parms.getImage("image");
    if (pixels == null) pixels = new Image3f(5,5,one);
    width  = (unsigned) pixels->width;
    height = (unsigned) pixels->height;

    Array2D<float> importance(height,width);
    for (size_t y = 0; y < height; y++)
      for (size_t x = 0; x < width; x++)
        importance.set(y, x, sinf(float(pi)*(y+0.5f)*rcp(float(height))) * reduce_add(pixels->get(x,y)));

    distribution = new Distribution2D(importance,width,height);
  }

  __forceinline Color HDRILight::Le(const Vector3f& wo) const
  {
    Vector3f wi = xfmVector(world2local, -wo);
    float theta = acosf(clamp(wi.y,-1.0f,1.0f));
    float phi = atan2f(-wi.z,-wi.x);
    if (phi < 0) phi += 2.0f * float(pi);
    float u = 1.0f - (phi * float(one_over_two_pi));
    float v = theta * float(one_over_pi);

    ssize_t x = clamp(ssize_t(u*width), ssize_t(0), ssize_t(width-1));
    ssize_t xNext = x + 1;
    if (size_t(xNext) == width) xNext = 0;
    float alpha = u*width - x;

    ssize_t y = clamp(ssize_t(v*height), ssize_t(0), ssize_t(height-1));
    ssize_t yNext = y+1;
    if (size_t(yNext) == height) yNext = height-1;
    float beta  = v*height - y;

    Color c0 = pixels->get(x,     y    );
    Color c1 = pixels->get(xNext, y    );
    Color c2 = pixels->get(xNext, yNext);
    Color c3 = pixels->get(x,     yNext);

    Color temp0 = beta*c3 + (1-beta)*c0;
    Color temp1 = beta*c2 + (1-beta)*c1;

    return L * (alpha*temp1 + (1-alpha)*temp0);
  }

  Color HDRILight::eval(const DifferentialGeometry& dg, const Vector3f& wi) const {
    return Le(-wi);
  }

  Color HDRILight::sample(const DifferentialGeometry& dg, Sample3f& wi, float& tMax, const Vec2f& sample) const
  {
    Sample2f pixelF = distribution->sample(sample);
    float theta = float(pi) * pixelF.value.y*rcp(float(height));
    float phi   = float(two_pi) * (1.0f - pixelF.value.x*rcp(float(width)));
    Vector3f _wi = Vector3f(-sinf(theta)*cosf(phi),cosf(theta),-sinf(theta)*sinf(phi));
    wi = Sample3f(xfmVector(local2world, _wi),pixelF.pdf*rcp(float(two_pi) * float(pi) * sinf(theta)));
    tMax = inf;
    return L*pixels->get(clamp(ssize_t(pixelF.value.x), ssize_t(0), ssize_t(width-1)),
                         clamp(ssize_t(pixelF.value.y), ssize_t(0), ssize_t(height-1)));
  }

  float HDRILight::pdf(const DifferentialGeometry& dg, const Vector3f& _wi) const {
    Vector3f wi = xfmVector(world2local, _wi);
    float theta = acosf(wi.y);
    float phi = atan2f(-wi.z,-wi.x);
    if (phi < 0) phi += 2.0f * float(pi);
    float u = 1.0f - float(one_over_two_pi) * phi;
    float v = float(one_over_pi) * theta;
    return distribution->pdf(Vec2f(u,v))*rcp(float(two_pi) * float(pi) * sinf(theta));
  }
}

