// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#ifndef __EMBREE_SWAPCHAIN_H__
#define __EMBREE_SWAPCHAIN_H__

#include <vector>
#include "framebuffer.h"

namespace embree
{
  /*! A swapchain is a sequence of framebuffers and a state. */
  template<typename BufferState = int> 
  class SwapChain : public RefCount
  {
  public:
    
    /*! constructs a swapchain of specified width, height, and number of buffers */
    SwapChain (size_t width, size_t height, size_t depth, const BufferState& initialState = BufferState())
      : width(width), height(height), depth(depth), buf(0)
    {
      for (size_t i=0; i<depth; i++) {
        _state.push_back(initialState);
        _buffer.push_back(new FrameBuffer(width,height));
      }
      _accu = new FrameBuffer(width,height);
    }
    
    /*! return the width of the swapchain */
    __forceinline size_t getWidth() const { return width;  }
    
    /*! return the height of the swapchain */
    __forceinline size_t getHeight() const { return height; }
    
    /*! return the number of buffers of the swapchain */
    __forceinline size_t getDepth() const { return depth;  }
    
    /*! goto next buffers */
    void swapBuffers() {
      buf = (buf+1)%depth;
      _buffer[buf]->wait();
    }

    /*! accumulate current framebuffer with accumulation buffer */
    void accumulate (Ref<FrameBuffer>& framebuffer, Vec2i start, Vec2i end, bool accumulate) 
    {
      if (accumulate) 
        for (ssize_t y=start.y; y<=end.y; y++)
          for (ssize_t x=start.x; x<=end.x; x++)
            framebuffer->get(x,y) = accu()->get(x,y) += framebuffer->get(x,y);
      else
        for (ssize_t y=start.y; y<=end.y; y++)
          for (ssize_t x=start.x; x<=end.x; x++)
            framebuffer->get(x,y) = accu()->get(x,y) = framebuffer->get(x,y);
    }

    /*! accumulate current framebuffer with accumulation buffer */
    void accumulateAndNormalize (Ref<FrameBuffer>& framebuffer, Vec2i start, Vec2i end, bool accumulate) 
    {
      if (accumulate) 
        for (ssize_t y=start.y; y<=end.y; y++) {
          for (ssize_t x=start.x; x<=end.x; x++) {
            const Vec4f a = accu()->get(x,y) += framebuffer->get(x,y);
            const Col3f c = Col3f(a.x,a.y,a.z) * rcp(a.w);
            framebuffer->get(x,y) = Vec4f(c.r,c.g,c.b,1.0f);
          }
        }
      else
        for (ssize_t y=start.y; y<=end.y; y++) {
          for (ssize_t x=start.x; x<=end.x; x++) {
            const Vec4f a = accu()->get(x,y) = framebuffer->get(x,y);
            const Col3f c = Col3f(a.x,a.y,a.z) * rcp(a.w);
            framebuffer->get(x,y) = Vec4f(c.r,c.g,c.b,1.0f);
          }
        }
    }
  
    /*! returns ID of current buffer */
    __forceinline size_t id() const { return buf; }

    /*! returns the accumulation buffer */
    __forceinline Ref<FrameBuffer>& accu() { return _accu; }
    
    /*! returns the current framebuffer */
    __forceinline Ref<FrameBuffer>& buffer() { return _buffer[buf]; }

    /*! returns the specified framebuffer */
    __forceinline Ref<FrameBuffer>& buffer(size_t id) { return _buffer[id]; }

    /*! returns the current state */
    __forceinline const BufferState& state() const { return _state[buf]; }
    __forceinline       BufferState& state()       { return _state[buf]; }

    /*! returns the specified state */
    __forceinline const BufferState& state(size_t id) const { return _state[id]; }
    __forceinline       BufferState& state(size_t id)       { return _state[id]; }
    
  private:
    size_t width;                             //!< width of the swapchain in pixels
    size_t height;                            //!< height of the swapchain in pixels
    size_t depth;                             //!< number of buffers in the swapchain
    size_t buf;                               //!< next buffer
    
  private:
    Ref<FrameBuffer> _accu;                   //!< special accumulation buffer
    std::vector<BufferState> _state;          //!< additional state to pass along the swapchain
    std::vector<Ref<FrameBuffer> > _buffer;   //!< the swapchain frame buffers
  };
}

#endif
