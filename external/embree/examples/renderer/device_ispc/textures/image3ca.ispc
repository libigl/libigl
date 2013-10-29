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

#include "image.isph"

struct Image3ca
{
  Image base;
  uint* data;
  bool owned;
};

inline uniform vec3f Image3ca__get_nearest_uniform(const uniform Image *uniform _this, const uniform int x, const uniform int y) 
{ 
  const uniform Image3ca *uniform this = (const uniform Image3ca *uniform)_this;
  const uniform float one_over_255 = (1.f/255.f);
  const uniform uint c = this->data[x+this->base.size.x*y]; 
  const uniform uint r = c         & 255;
  const uniform uint g = (c >>  8) & 255;
  const uniform uint b = (c >> 16) & 255;
  return make_vec3f(r*one_over_255, g*one_over_255, b*one_over_255);
}

inline varying vec3f Image3ca__get_nearest_varying(const uniform Image *uniform _this, const varying int x, const varying int y) 
{ 
  const uniform Image3ca *uniform this = (const uniform Image3ca *uniform)_this;
  const uniform float one_over_255 = (1.f/255.f);
  const varying uint idx = x+this->base.size.x*y;
  const varying uint c = this->data[idx]; 
  const varying uint r = c         & 255;
  const varying uint g = (c >>  8) & 255;
  const varying uint b = (c >> 16) & 255;
  return make_vec3f(r*one_over_255, g*one_over_255, b*one_over_255);
}

varying vec3f Image3ca__get_bilinear_varying(const uniform Image *uniform _this, const varying float x, const varying float y) 
{ 
  const uniform Image3ca *uniform this = (const uniform Image3ca *uniform)_this;

  /* repeat texture */
  const float xr = (x-floor(x))*(float)this->base.size.x;
  const float yr = (y-floor(y))*(float)this->base.size.y;

  /* calculate interpolation weights */
  const float sx = xr-floor(xr);
  const float sy = yr-floor(yr);

  /* read texels */
  const int x0 = clamp((int)xr,0,this->base.size.x-1);
  const int y0 = clamp((int)yr,0,this->base.size.y-1);
  const int x1 = x0+1 == this->base.size.x ? 0 : x0+1;
  const int y1 = y0+1 == this->base.size.y ? 0 : y0+1;

  const vec3f c00 = Image3ca__get_nearest_varying(_this,x0,y0);
  const vec3f c01 = Image3ca__get_nearest_varying(_this,x0,y1);
  const vec3f c10 = Image3ca__get_nearest_varying(_this,x1,y0);
  const vec3f c11 = Image3ca__get_nearest_varying(_this,x1,y1);
  
  /* interpolate texels */
  return lerp(sy,
              lerp(sx,c00,c10),
              lerp(sx,c01,c11));
}

void Image3ca__Destructor(uniform RefCount* uniform _this)
{
  uniform Image3ca* uniform this = (uniform Image3ca* uniform) _this;
  if (this->owned) delete[] this->data;
  Image__Destructor(_this);
}

void Image3ca__Constructor(uniform Image3ca* uniform this, const uniform vec2ui size, 
						   uniform uint *uniform data, uniform bool owned)
{
  Image__Constructor(&this->base,
                     Image3ca__Destructor,size,
                     Image3ca__get_nearest_uniform,
                     Image3ca__get_nearest_varying,
                     Image3ca__get_bilinear_varying);
  this->data = data;
  this->owned = owned;
}

export void* uniform Image3ca__new(const uniform int size_x, const uniform int size_y, uniform uint *uniform data, const uniform bool copy)
{
  uniform Image3ca *uniform this = uniform new uniform Image3ca;
  if (copy) {
    uniform uint* uniform local = uniform new uniform uint[size_x*size_y];
    memcpy(local,data,size_x*size_y*sizeof(uniform uint));
    Image3ca__Constructor(this,make_vec2ui(size_x,size_y),local,copy);
  } else {
    Image3ca__Constructor(this,make_vec2ui(size_x,size_y),data,copy);
  }
  return this;
}
