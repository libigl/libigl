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

#include "swapchain.isph"
#include "framebuffer_rgb_float32.isph"
#include "framebuffer_rgba8.isph"
#include "framebuffer_rgb8.isph"

typedef uniform FrameBuffer* uniform FrameBuffer_uptr;
typedef uniform FrameBuffer* uniform (*FrameBuffer__new)(const uniform uint width, const uniform uint height, const void* uniform ptr);

void SwapChain__Destructor(uniform RefCount* uniform _this)
{
  LOG(print("SwapChain__Destructor\n"));
  uniform SwapChain* uniform this = (uniform SwapChain* uniform) _this;

  for (uniform int i=0; i<this->depth; i++) {
    RefCount__Destroy(&this->buffers[i]->base); 
    this->buffers[i] = NULL;
  }
  delete[] this->buffers; this->buffers = NULL;

  RefCount__Destroy(&this->accu->base);
  RefCount__Destructor(&this->base);
}

void SwapChain__Constructor(uniform SwapChain* uniform this,
                            const uniform uint width, 
                            const uniform uint height,
                            const uniform uint depth,
                            const void* uniform* uniform ptrs,
                            const uniform DestructorFunc destructor,
                            const uniform FrameBuffer__new create)
{
  LOG(print("SwapChain__Constructor\n"));
  RefCount__Constructor(&this->base,destructor);
  this->width = width;
  this->height = height;
  this->depth = depth;
  this->buf = 0;
  this->buffers = uniform new FrameBuffer_uptr [depth];
  for (uniform int i=0; i<depth; i++) {
    if (ptrs) this->buffers[i] = create(width,height,ptrs[i]);
    else      this->buffers[i] = create(width,height,NULL);
  }
  this->accu = AccuBuffer__new(width,height);
}

export void* uniform SwapChainRGBFloat32__new(const uniform uint width, const uniform uint height, const uniform uint depth, const void* uniform* uniform ptrs)
{
  uniform SwapChain* uniform this = uniform new uniform SwapChain;
  SwapChain__Constructor(this,width,height,depth,ptrs,SwapChain__Destructor,FrameBufferRGBFloat32__new);
  return this;
}

export void* uniform SwapChainRGBA8__new(const uniform uint width, const uniform uint height, const uniform uint depth, const void* uniform* uniform ptrs)
{
  uniform SwapChain* uniform this = uniform new uniform SwapChain;
  SwapChain__Constructor(this,width,height,depth,ptrs,SwapChain__Destructor,FrameBufferRGBA8__new);
  return this;
}

#if !defined(__MIC__)
export void* uniform SwapChainRGB8__new(const uniform uint width, const uniform uint height, const uniform uint depth, const void* uniform* uniform ptrs)
{
  uniform SwapChain* uniform this = uniform new uniform SwapChain;
  SwapChain__Constructor(this,width,height,depth,ptrs,SwapChain__Destructor,FrameBufferRGB8__new);
  return this;
}
#endif

uniform FrameBuffer* uniform SwapChain__get_buffer(uniform SwapChain* uniform this) {
  return this->buffers[this->buf];
}

uniform AccuBuffer* uniform SwapChain__get_accu(uniform SwapChain* uniform this) {
  return this->accu;
}

export void* uniform SwapChain__map(void* uniform _this, uniform int id)
{
  uniform SwapChain* uniform this = (uniform SwapChain* uniform) _this;
  if (id < 0) id = this->buf;
  uniform FrameBuffer* uniform fb = this->buffers[id];
  return fb->map(fb);
}

export void SwapChain__unmap(void* uniform _this)
{
}

export void SwapChain__swap(void* uniform _this)
{
  uniform SwapChain* uniform this = (uniform SwapChain* uniform) _this;
  this->buf = (this->buf+1) % this->depth;
}
