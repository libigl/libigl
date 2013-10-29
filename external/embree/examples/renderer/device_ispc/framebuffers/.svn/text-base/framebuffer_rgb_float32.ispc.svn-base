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

#include "framebuffer_rgb_float32.isph"

void FrameBufferRGBFloat32__set(uniform FrameBuffer* uniform _this, const int x, const int y, const vec3f rgb) 
{ 
  uniform FrameBufferRGBFloat32* uniform this = (uniform FrameBufferRGBFloat32* uniform) _this;
  this->ptr[x+this->base.size.x*y] = rgb; 
}

void* uniform FrameBufferRGBFloat32__map(uniform FrameBuffer* uniform _this) 
{
  uniform FrameBufferRGBFloat32* uniform this = (uniform FrameBufferRGBFloat32* uniform) _this;
  return this->ptr;
}

void FrameBufferRGBFloat32__Destructor(uniform RefCount* uniform _this)
{
  uniform FrameBufferRGBFloat32* uniform this = (uniform FrameBufferRGBFloat32* uniform) _this;
  if (this->allocated) delete[] this->ptr; this->ptr = NULL;
  FrameBuffer__Destructor(&this->base);
}

void FrameBufferRGBFloat32__Constructor(uniform FrameBufferRGBFloat32* uniform this, const uniform uint width, const uniform uint height, const void* uniform ptr)
{
  FrameBuffer__Constructor(&this->base,width,height,
                           FrameBufferRGBFloat32__Destructor,
                           FrameBufferRGBFloat32__set,
                           FrameBufferRGBFloat32__map);
  if (ptr) {
    this->allocated = false;
    this->ptr = (uniform vec3f* uniform) ptr;
  } else {
    this->allocated = true;
    this->ptr = uniform new uniform vec3f[width*height];
  }
  memset(this->ptr,0,width*height*sizeof(uniform vec3f));
}

uniform FrameBuffer* uniform FrameBufferRGBFloat32__new(const uniform uint width, const uniform uint height, const void* uniform ptr)
{
  uniform FrameBufferRGBFloat32* uniform this = uniform new uniform FrameBufferRGBFloat32;
  FrameBufferRGBFloat32__Constructor(this,width,height,ptr);
  return (uniform FrameBuffer* uniform) this;
}
