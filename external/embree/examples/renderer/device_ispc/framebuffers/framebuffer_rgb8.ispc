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

#include "framebuffer_rgb8.isph"

void FrameBufferRGB8__set(uniform FrameBuffer* uniform _this, const int x, const int y, const vec3f rgb) 
{ 
  uniform FrameBufferRGB8* uniform this = (uniform FrameBufferRGB8* uniform) _this;
  int idx = x+this->base.size.x*y;
  this->ptr[idx].x = (int) clamp(rgb.x * 255.0f, 0.0f, 255.0f);
  this->ptr[idx].y = (int) clamp(rgb.y * 255.0f, 0.0f, 255.0f);
  this->ptr[idx].z = (int) clamp(rgb.z * 255.0f, 0.0f, 255.0f);
}

void* uniform FrameBufferRGB8__map(uniform FrameBuffer* uniform _this) 
{
  uniform FrameBufferRGB8* uniform this = (uniform FrameBufferRGB8* uniform) _this;
  return this->ptr;
}

void FrameBufferRGB8__Destructor(uniform RefCount* uniform _this)
{
  uniform FrameBufferRGB8* uniform this = (uniform FrameBufferRGB8* uniform) _this;
  if (this->allocated) delete[] this->ptr; this->ptr = NULL;
  FrameBuffer__Destructor(&this->base);
}

void FrameBufferRGB8__Constructor(uniform FrameBufferRGB8* uniform this, const uniform uint width, const uniform uint height, const void* uniform ptr)
{
  FrameBuffer__Constructor(&this->base,width,height,
                           FrameBufferRGB8__Destructor,
                           FrameBufferRGB8__set,
                           FrameBufferRGB8__map);
  if (ptr) {
    this->allocated = false;
    this->ptr = (uniform vec3uc* uniform) ptr;
  } else {
    this->allocated = true;
    this->ptr = uniform new uniform vec3uc[width*height];
  }
  memset(this->ptr,0,width*height*sizeof(uniform vec3uc));
}

uniform FrameBuffer* uniform FrameBufferRGB8__new(const uniform uint width, const uniform uint height, const void* uniform ptr)
{
  uniform FrameBufferRGB8* uniform this = uniform new uniform FrameBufferRGB8;
  FrameBufferRGB8__Constructor(this,width,height,ptr);
  return (uniform FrameBuffer* uniform) this;
}
