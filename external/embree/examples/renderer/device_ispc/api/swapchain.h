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

#ifndef __EMBREE_SWAPCHAIN_H__
#define __EMBREE_SWAPCHAIN_H__

#include <vector>
#include "framebuffer.h"

namespace embree
{
  extern "C" int g_serverCount;
  extern "C" int g_serverID;

  /*! A swapchain is a sequence of framebuffers and a state. */
  class SwapChain : public RefCount
  {
  public:
    typedef FrameBuffer* (*Factory)(size_t width, size_t height);

  public:
    
    /*! constructs a swapchain of specified width, height, and number of buffers */
    SwapChain (const std::string format, size_t width, size_t height, size_t depth, Factory fb_factory)
      : format(format), width(width), height(height), depth(depth), buf(0)
    {
      for (size_t i=0; i<depth; i++) {
        _buffer.push_back(fb_factory(width,height));
      }
      _accu = new AccuBuffer(width,height);
    }

    /*! return the width of the swapchain */
    __forceinline size_t getWidth() const { return width;  }
    
    /*! return the height of the swapchain */
    __forceinline size_t getHeight() const { return height; }
    
    /*! return the number of buffers of the swapchain */
    __forceinline size_t getDepth() const { return depth;  }

    /*! determines if a line is active when rendering in network mode */
    __forceinline bool activeLine(int py) {
      const int row = py>>2;
      return ((row+g_serverID) % g_serverCount) == 0;
    }
    
    /*! determines the first active line */
    __forceinline size_t firstActiveLine() {
      return g_serverID;
    }

    /*! determines the destination line in network mode */
    __forceinline size_t raster2buffer(size_t py) {
      return 4*((py>>2)/g_serverCount)+(py&0x3);
    }
    
    /*! goto next buffer */
    void swapBuffers() {
      buf = (buf+1)%depth;
      _buffer[buf]->wait();
    }

    /*! Clear a part of the framebuffer. */
    void clear(Vec2i start, Vec2i end) 
    {
      for (ssize_t y=start.y; y<=end.y; y++) {
        for (ssize_t x=start.x; x<=end.x; x++) {
          _accu->clear(x,y);
        }
      }
    }

    /*! accumulate inside accumulation buffer */
    void accumulate(size_t x, size_t y, const Color& color, const float weight) {
      _accu->add(x,y,Vec4f(color.r,color.g,color.b,weight));
    }

    /*! accumulate inside accumulation buffer */
    void update(size_t x, size_t y, const Color& color, const float weight, const bool accumulate) {
      _accu->update(x,y,Vec4f(color.r,color.g,color.b,weight),accumulate);
    }

    /*! returns framebuffer format */
    __forceinline std::string getFormat() const { return format; }

    /*! returns accumulation buffer */
    __forceinline Ref<AccuBuffer>& accu() { return _accu; }

    /*! returns ID of current buffer */
    __forceinline size_t id() const { return buf; }

    /*! returns the current framebuffer */
    __forceinline Ref<FrameBuffer>& buffer() { return _buffer[buf]; }

    /*! returns the specified framebuffer */
    __forceinline Ref<FrameBuffer>& buffer(size_t id) { return _buffer[id]; }

  private:
    std::string format;
    size_t width;                             //!< width of the swapchain in pixels
    size_t height;                            //!< height of the swapchain in pixels
    size_t depth;                             //!< number of buffers in the swapchain
    size_t buf;                               //!< next buffer
    
  private:
    Ref<AccuBuffer> _accu;                    //!< special accumulation buffer
    std::vector<Ref<FrameBuffer> > _buffer;   //!< the swapchain frame buffers
  };
}

#endif
