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

#include "texture.isph"
#include "textures/image.isph"

struct NearestNeighbor
{
  Texture base;
  Image* image;
};

varying vec3f NearestNeighbor__get(const uniform Texture *uniform _this, const varying vec2f& p)
{
  const uniform NearestNeighbor* uniform this = (const uniform NearestNeighbor* uniform) _this;
  const uniform Image* uniform image = this->image;
  const uniform vec2ui size = image->size;
  const float x1 = p.x-floor(p.x);
  const float y1 = p.y-floor(p.y);
  const float x2 = x1*(float)size.x;
  const float y2 = y1*(float)size.y;
  const int x3 = clamp((int)x2,0,size.x-1);
  const int y3 = clamp((int)y2,0,size.y-1);
  return image->get_nearest_varying(image,x3,y3);
}

void NearestNeighbor__Destructor(uniform RefCount* uniform _this)
{
  uniform NearestNeighbor* uniform this = (uniform NearestNeighbor* uniform) _this;
  RefCount__DecRef(&this->image->base);
  Texture__Destructor(_this);
}

void NearestNeighbor__Constructor(uniform NearestNeighbor *uniform this, uniform Image* uniform image)
{
  Texture__Constructor(&this->base,NearestNeighbor__Destructor,NearestNeighbor__get);
  RefCount__IncRef(&image->base);
  this->image = image;
}

export void* uniform NearestNeighbor__new(void* uniform image)
{
  uniform NearestNeighbor* uniform tex = uniform new uniform NearestNeighbor;
  NearestNeighbor__Constructor(tex,(uniform Image* uniform)image);
  return tex;
}
