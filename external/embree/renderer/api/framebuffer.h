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

#ifndef __EMBREE_FRAMEBUFFER_H__
#define __EMBREE_FRAMEBUFFER_H__

#include "default.h"

namespace embree
{
  /*! Normalizes a color. */
  __forceinline Vec4f normalizeColor(const Vec4f& c) { 
    return c * rcp(c.w); 
  }
  
  /*! framebuffers the swapchain consists of */
  struct FrameBuffer : public RefCount
  {
    /*! constructs a new framebuffer of specified size */
    FrameBuffer (size_t width, size_t height)
      : width(width), height(height)
    {
      data = new Vec4f[width*height];
      memset(data,0,width*height*sizeof(Vec4f));
      remainingTiles = 0;
    }
    
    /*! destroys the framebuffer */
    ~FrameBuffer () {
      delete[] data; data = NULL;
    }
    
    /*! return the width of the swapchain */
    __forceinline size_t getWidth() const { return width;  }
    
    /*! return the height of the swapchain */
    __forceinline size_t getHeight() const { return height; }

    /*! return a pointer to the raw pixel data */
    __forceinline float* getData() { return (float*)data; }

    /*! signals the framebuffer that rendering starts */
    void startRendering(size_t numTiles = 1) {
      Lock<MutexSys> lock(mutex);
      remainingTiles = numTiles;
    }
    
    /*! register a tile as finished */
    bool finishTile()
    {
      Lock<MutexSys> lock(mutex);
      remainingTiles--;
      if (remainingTiles == 0) {
        condition.broadcast();
        return true;
      }
      return false;
    }
    
    /*! wait for rendering to finish */
    void wait() {
      Lock<MutexSys> lock(mutex);
      while (remainingTiles != 0) condition.wait(mutex);
    }
    
    /*! access to pixels */
    __forceinline const Vec4f& get(size_t x, size_t y) const { return data[y*width+x]; }
    __forceinline       Vec4f& get(size_t x, size_t y)       { return data[y*width+x]; }
    
    /*! Clear a part of the framebuffer. */
    void clear(Vec2i start, Vec2i end) {
      for (ssize_t y=start.y; y<=end.y; y++) {
        for (ssize_t x=start.x; x<=end.x; x++) {
          data[y*width+x] = Vec4f(0.0f,0.0f,0.0f,1E-10f);
        }
      }
    }
    
    /*! Add a new sample to the framebuffer. \param p is the location of the
     *  sample \param start is the beginning of the tile \param end is
     *  the end of the tile \param color is the color of the sample
     *  \param weight is the weight of the sample */
    void accumulate(const Vec2f& p, Vec2i start, Vec2i end, const Col3f& color, const float weight)
    {
      /*! ignore samples outside the tile */
      int x = (int)p.x, y = (int)p.y;
      if (x < start.x || x > end.x || y < start.y || y > end.y) return;
      
      /*! accumulate color and weight */
      data[y*width+x] += Vec4f(color.r,color.g,color.b,weight);
    }
    
  private:
    size_t width;              //!< width of the framebuffer in pixels
    size_t height;             //!< height of the framebuffer in pixels
    Vec4f* data;               //!< framebuffer data

  private:
    size_t remainingTiles;     //!< number of tiles that are not rendered yet
    MutexSys mutex;            //!< mutex to protect access to remainingTiles
    ConditionSys condition;    //!< condition to signal threads waiting for render to finish
  };
}

#endif
