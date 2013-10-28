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

#ifndef __EMBREE_FRAMEBUFFER_H__
#define __EMBREE_FRAMEBUFFER_H__

#include "../default.h"

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
    FrameBuffer (size_t width, size_t height, size_t stride, void* ptr)
      : width(width), height(height), stride(stride), data(ptr), allocated(false), remainingTiles(0) {}

    /*! read pixel */
    virtual const Color get(size_t x, size_t y) const = 0;

    /*! write pixel */
    virtual void set(size_t x, size_t y, const Color& c) = 0;
    
    /*! return the width of the framebuffer */
    __forceinline size_t getWidth() const { return width;  }
    
    /*! return the height of the framebuffer */
    __forceinline size_t getHeight() const { return height; }

    /*! return the stride of the framebuffer */
    __forceinline size_t getStride() const { return stride; }

    /*! return a pointer to the raw pixel data */
    __forceinline void* getData() { return data; }

    /*! signals the framebuffer that rendering starts */
    void startRendering(size_t numTiles = 1) {
      Lock<MutexSys> lock(mutex);
      remainingTiles = numTiles;
    }
    
    /*! register a tile as finished */
    bool finishTile(int numTiles = 1)
    {
      Lock<MutexSys> lock(mutex);
      remainingTiles-=numTiles;
      if (remainingTiles == 0) {
        condition.broadcast();
        return true;
      }
      return false;
    }
    
    /*! wait for rendering to finish */
    virtual void wait() {
      Lock<MutexSys> lock(mutex);
      while (remainingTiles != 0) condition.wait(mutex);
    }
    
  protected:
    size_t width;              //!< width of the framebuffer in pixels
    size_t height;             //!< height of the framebuffer in pixels
    size_t stride;             //!< stride of line on bytes
    void* data;                //!< framebuffer data
    bool allocated;            //!< true if framebuffer allocated by us

  private:
    size_t remainingTiles;     //!< number of tiles that are not rendered yet
    MutexSys mutex;            //!< mutex to protect access to remainingTiles
    ConditionSys condition;    //!< condition to signal threads waiting for render to finish
  };

  /*! RGB_FLOAT32 framebuffer */
  struct FrameBufferRGBFloat32 : public FrameBuffer
  {
  public:

    /*! class factory */
    static FrameBuffer* create(size_t width, size_t height, void* ptr) {
      return new FrameBufferRGBFloat32(width,height,ptr);
    }

  protected:

    /*! constructs a new framebuffer of specified size */
    FrameBufferRGBFloat32 (size_t width, size_t height, void* ptr) 
      : FrameBuffer(width,height,width*sizeof(Col3f),ptr) 
    {
      if (!data) {
        data = (void*) new Col3f[width*height];
        allocated = true;
      }
      memset(data,0,width*height*sizeof(Col3f));
    }
    
    /*! destroys the framebuffer */
    ~FrameBufferRGBFloat32 () {
      if (allocated) delete[] (Col3f*) data; data = NULL;
    }

    /*! read pixel */
    const Color get(size_t x, size_t y) const {
      Col3f c = ((Col3f*)data)[y*width+x];
      return Color(c.r,c.g,c.b);
    }

    /*! write pixel */
    void set(size_t x, size_t y, const Color& c) {
      ((Col3f*)data)[y*width+x] = Col3f(c.r,c.g,c.b);
    }
  };

  /*! RGBA8 framebuffer */
  struct FrameBufferRGBA8 : public FrameBuffer
  {
  public:

    /*! class factory */
    static FrameBuffer* create(size_t width, size_t height, void* ptr) {
      return new FrameBufferRGBA8(width,height,ptr);
    }

  protected:

    /*! constructs a new framebuffer of specified size */
    FrameBufferRGBA8 (size_t width, size_t height, void* ptr)
      : FrameBuffer(width,height,4*width,ptr) 
    {
      if (!data) {
        data = malloc(stride*height);
        allocated = true;
      }
      memset(data,0,stride*height);
    }
    
    /*! destroys the framebuffer */
    ~FrameBufferRGBA8 () {
      if (allocated) free(data); data = NULL;
    }

    /*! read pixel */
    const Color get(size_t x, size_t y) const 
    {
      const float one_over_255 = 1.0f/255.0f;
      unsigned char* pixel = (unsigned char*)data+stride*y+4*x;
      return Color(pixel[0]*one_over_255,
                   pixel[1]*one_over_255,
                   pixel[2]*one_over_255);
    }

    /*! write pixel */
    void set(size_t x, size_t y, const Color& c) 
    {
      unsigned char* pixel = (unsigned char*)data+stride*y+4*x;
      pixel[0] = (unsigned char) clamp(c.r*255.0f,0.0f,255.0f);
      pixel[1] = (unsigned char) clamp(c.g*255.0f,0.0f,255.0f);
      pixel[2] = (unsigned char) clamp(c.b*255.0f,0.0f,255.0f);
      pixel[3] = 0;
    }
  };

  /*! RGB8 framebuffer */
  struct FrameBufferRGB8 : public FrameBuffer
  {
  public:

    /*! class factory */
    static FrameBuffer* create(size_t width, size_t height, void* ptr) {
      return new FrameBufferRGB8(width,height,ptr);
    }

  protected:

    /*! constructs a new framebuffer of specified size */
    FrameBufferRGB8 (size_t width, size_t height, void* ptr)
      : FrameBuffer(width,height,(3*width+3)/4*4,ptr) 
    {
      if (!data) {
        data = malloc(stride*height);
        allocated = true;
      }
      memset(data,0,stride*height);
    }
    
    /*! destroys the framebuffer */
    ~FrameBufferRGB8 () {
      if (allocated) free(data); data = NULL;
    }

    /*! read pixel */
    const Color get(size_t x, size_t y) const 
    {
      const float one_over_255 = 1.0f/255.0f;
      unsigned char* pixel = (unsigned char*)data+stride*y+3*x;
      return Color(pixel[0]*one_over_255,
                   pixel[1]*one_over_255,
                   pixel[2]*one_over_255);
    }

    /*! write pixel */
    void set(size_t x, size_t y, const Color& c) 
    {
      unsigned char* pixel = (unsigned char*)data+stride*y+3*x;
      pixel[0] = (unsigned char) clamp(c.r*255.0f,0.0f,255.0f);
      pixel[1] = (unsigned char) clamp(c.g*255.0f,0.0f,255.0f);
      pixel[2] = (unsigned char) clamp(c.b*255.0f,0.0f,255.0f);
    }
  };

  /*! accumulation buffer */
  struct AccuBuffer : public RefCount
  {
  public:

    /*! constructs a new framebuffer of specified size */
    AccuBuffer (size_t width, size_t height) 
    : width(width), height(height), data(NULL)
    {
      data = new Vec4f[width*height];
      memset(data,0,width*height*sizeof(Vec4f));
    }
    
    /*! destroys the framebuffer */
    ~AccuBuffer () {
      delete[] data; data = NULL;
    }

    /*! return the width of the swapchain */
    __forceinline size_t getWidth() const { return width;  }
    
    /*! return the height of the swapchain */
    __forceinline size_t getHeight() const { return height; }

    /*! return a pointer to the raw pixel data */
    __forceinline void* getData() { return data; }

    /*! clear buffer */
    __forceinline void clear(size_t x, size_t y) {
      data[y*width+x] = Vec4f(0.0f,0.0f,0.0f,1E-10f);
    }

    /*! set pixel */
    __forceinline void set(size_t x, size_t y, const Vec4f& c) {
      data[y*width+x] = c;
    }

    /*! accumulate pixel */
    __forceinline void add(size_t x, size_t y, const Vec4f& c) {
      data[y*width+x] += c;
    }

    /*! update pixel */
    __forceinline Color update(size_t x, size_t y, const Color& c, const float weight, bool accu) 
    {
      if (accu) {
        const Vec4f cur = data[y*width+x];
        const Vec4f next = cur + Vec4f(c.r,c.g,c.b,weight);
        data[y*width+x] = next;
        const float norm = rcp(next.w);
        return Color(next.x,next.y,next.z)*norm;
      }
      else {
        data[y*width+x] = Vec4f(c.r,c.g,c.b,weight);
        return c*rcp(weight);
      }
    }

    /*! read pixel */
    __forceinline const Color get(size_t x, size_t y) const 
    {
      const Vec4f& c = data[y*width+x];
      const float norm = rcp(c.w);
      return Color(c.x,c.y,c.z)*norm;
    }

  protected:
    size_t width;              //!< width of the framebuffer in pixels
    size_t height;             //!< height of the framebuffer in pixels
    Vec4f* data;                //!< framebuffer data
  };
}

#endif
