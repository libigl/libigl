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

#ifndef __EMBREE_IMAGE_H__
#define __EMBREE_IMAGE_H__

#include <vector>
#include <string>

#include "sys/platform.h"
#include "sys/ref.h"
#include "sys/filename.h"
#include "math/col3.h"

namespace embree
{
  /* virtual interface to image */
  class Image : public RefCount {
  public:
    Image (size_t width, size_t height, const std::string& name) : width(width), height(height), name(name) {}
    virtual ~Image() {}
    virtual Col3f get(size_t x, size_t y) const = 0;
    virtual void  set(size_t x, size_t y, const Col3f& c) = 0;
  public:
    size_t width,height;
    std::string name;
  };

  /* conversion to and from Col3f */
  __forceinline Col3f toCol3f(const Col3f& c) { return c; }
  __forceinline Col3f toCol3f(const Col3c& c) { return Col3f(c.r,c.g,c.b)*(1.0f/255.0f); }

  __forceinline void setCol3f(Col3f& d, const Col3f& c) { d = c; }
  __forceinline void setCol3f(Col3c& d, const Col3f& c) { d = Col3c(char(clamp(c.r)*255.0f),char(clamp(c.g)*255.0f),char(clamp(c.b)*255.0f)); }

  /* main image class templated over element type */
  template<typename T>
    class ImageT : public Image {
  public:
    
    /*! create empty image */
    ImageT (size_t width = 0, size_t height = 0, const std::string& name = "")
      : Image(width,height,name)
    {
      data = (T*) alignedMalloc(width*height*sizeof(T),64);
      for (size_t i=0; i<width*height; i++) data[i] = zero;
    }

    /*! create image of constant color */
    ImageT (size_t width, size_t height, const T& color, const std::string& name = "")
      : Image(width,height,name)
    {
      data = (T*) alignedMalloc(width*height*sizeof(T),64);
      for (size_t i=0; i<width*height; i++) data[i] = color;
    }

    /*! initialize image from color data */
    ImageT (size_t width, size_t height, const T* color, const std::string& name = "")
    : Image(width,height,name)
    {
      data = (T*) alignedMalloc(width*height*sizeof(T),64);
      for (size_t i=0; i<width*height; i++) data[i] = color[i];
    }

    /*! image destruction */
    virtual ~ImageT() {
      if (data) alignedFree(data); data = NULL;
    }
    
    /*! returns pixel color */
    __forceinline Col3f get(size_t x, size_t y) const { 
      return toCol3f(data[y*width+x]); 
    }

    /*! sets pixel */
    __forceinline void set(size_t x, size_t y, const Col3f& c) { 
      setCol3f(data[y*width+x],c); 
    }

    /*! returns data pointer of image */
    __forceinline void* ptr() {
      return (void*)data;
    }

  private:
    T* data;
  };
  
  /*! Shortcuts for common image types. */
  typedef ImageT<Col3c> Image3c;
  typedef ImageT<Col3f> Image3f;
  
  /*! Loads image from PPM file. */
  Ref<Image> loadPPM(const FileName& fileName);
  
  /*! Loads image from PFM file. */
  Ref<Image> loadPFM(const FileName& fileName);
  
  /*! Loads image from EXR file. */
  Ref<Image> loadExr(const FileName& fileName);
  
  /*! Loads image using ImageMagick. */
  Ref<Image> loadMagick(const FileName& fileName);
  
  /*! Loads image from file. Format is auto detected. */
  Ref<Image> loadImage(const FileName& filename, bool cache = false);
  
  /*! Store image to PPM file. */
  void storePPM(const Ref<Image>& img, const FileName& fileName);
  
  /*! Store image to PFM file. */
  void storePFM(const Ref<Image>& img, const FileName& fileName);
  
  /*! Store image to TGA file. */
  void storeTga(const Ref<Image>& img, const FileName& fileName);
  
  /*! Store image to EXR file. */
  void storeExr(const Ref<Image>& img, const FileName& fileName);
  
  /*! Store image to file using ImageMagick. */
  void storeMagick(const Ref<Image>& img, const FileName& fileName);
  
  /*! Store image to file. Format is auto detected. */
  void storeImage(const Ref<Image>& image, const FileName& filename);
}

#endif
