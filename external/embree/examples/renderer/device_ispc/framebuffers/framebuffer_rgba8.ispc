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

#include "framebuffer_rgba8.isph"

void FrameBufferRGBA8__set(uniform FrameBuffer* uniform _this, const int x, const int y, const vec3f rgb) 
{ 
  uniform FrameBufferRGBA8* uniform this = (uniform FrameBufferRGBA8* uniform) _this;
  const int index = x+this->base.size.x*y;
  int r = (int) clamp(rgb.x * 255.0f, 0.0f, 255.0f);
  int g = (int) clamp(rgb.y * 255.0f, 0.0f, 255.0f);
  int b = (int) clamp(rgb.z * 255.0f, 0.0f, 255.0f);
  this->ptr[index] = (b << 16) | (g << 8) | r; 
}

void* uniform FrameBufferRGBA8__map(uniform FrameBuffer* uniform _this) 
{
  uniform FrameBufferRGBA8* uniform this = (uniform FrameBufferRGBA8* uniform) _this;
  return this->ptr;
}

void FrameBufferRGBA8__Destructor(uniform RefCount* uniform _this)
{
  uniform FrameBufferRGBA8* uniform this = (uniform FrameBufferRGBA8* uniform) _this;
  if (this->allocated) delete[] this->ptr; this->ptr = NULL;
  FrameBuffer__Destructor(&this->base);
}

void FrameBufferRGBA8__Constructor(uniform FrameBufferRGBA8* uniform this, const uniform uint width, const uniform uint height, const void* uniform ptr)
{
  FrameBuffer__Constructor(&this->base,width,height,
                           FrameBufferRGBA8__Destructor,
                           FrameBufferRGBA8__set,
                           FrameBufferRGBA8__map);
  if (ptr) {
    this->allocated = false;
    this->ptr = (uniform int* uniform) ptr;
  } else {
    this->allocated = true;
    this->ptr = uniform new uniform int[width*height];
  }
  memset(this->ptr,0,width*height*sizeof(uniform int));
}

uniform FrameBuffer* uniform FrameBufferRGBA8__new(const uniform uint width, const uniform uint height, const void* uniform ptr)
{
  uniform FrameBufferRGBA8* uniform this = uniform new uniform FrameBufferRGBA8;
  FrameBufferRGBA8__Constructor(this,width,height,ptr);
  return (uniform FrameBuffer* uniform) this;
}
