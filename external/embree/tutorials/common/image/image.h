// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
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

#pragma once

#include <vector>
#include <string>

#include "sys/platform.h"
#include "sys/ref.h"
#include "sys/filename.h"
#include "math/color.h"

namespace embree
{
  /* virtual interface to image */
  class Image : public RefCount {
  public:
    Image (size_t width, size_t height, const std::string& name) : width(width), height(height), name(name) {}
    virtual ~Image() {}
    virtual Color4 get(size_t x, size_t y) const = 0;
    virtual void   set(size_t x, size_t y, const Color4& c) = 0;
    void set(size_t x, size_t y, const Color& c) { set(x,y,Color4(c.r,c.g,c.b,1.0f)); }
  public:
    size_t width,height;
    std::string name;
  };

  /* main image class templated over element type */
  template<typename T>
    class ImageT : public Image {
  public:
    
    /*! create empty image */
    ImageT (size_t width = 0, size_t height = 0, const std::string& name = "")
      : Image(width,height,name)
    {
      data = (T*) malloc(width*height*sizeof(T));
      memset(data,0,width*height*sizeof(T));
    }

    /*! create image of constant color */
    ImageT (size_t width, size_t height, const T& color, const std::string& name = "")
      : Image(width,height,name)
    {
      data = (T*) malloc(width*height*sizeof(T));
      for (size_t i=0; i<width*height; i++) data[i] = color;
    }

    /*! initialize image from color data */
    ImageT (size_t width, size_t height, T* color, const bool copy = true, const std::string& name = "")
      : Image(width,height,name)
    {
      if (copy) {
        data = (T*) malloc(width*height*sizeof(T));
        for (size_t i=0; i<width*height; i++) data[i] = color[i];
      } 
      else {
        data = color;
      } 
    }

    /*! image destruction */
    virtual ~ImageT() {
      if (data) free(data); data = NULL;
    }
    
    /*! returns pixel color */
    __forceinline Color4 get(size_t x, size_t y) const { 
      return Color4(data[y*width+x]); 
    }

    /*! sets pixel */
    __forceinline void set(size_t x, size_t y, const Color4& c) { 
      c.set(data[y*width+x]); 
    }

    /*! returns data pointer of image */
    __forceinline void* ptr() {
      return (void*)data;
    }

    /*! returns and forgets about data pointer of image */
    __forceinline void* steal_ptr() {
      T* ptr = data;
      data = NULL;
      return (void*)ptr;
    }

  protected:
    T* data;
  };
  
  /*! Shortcuts for common image types. */
  typedef ImageT<Col3c> Image3c;
  typedef ImageT<Col3f> Image3f;
  typedef ImageT<Col4c> Image4c;
  typedef ImageT<Col4f> Image4f;
  
  /*! Generate a JPEG encoded image from a RGB8 buffer in memory. */
  void encodeRGB8_to_JPEG(unsigned char *image, size_t width, size_t height, unsigned char **encoded, size_t *capacity);

  /*! Loads image from EXR file. */
  Ref<Image> loadExr(const FileName& fileName);
  
  /*! Loads image from file. Format is auto detected. */
  Ref<Image> loadImage(const FileName& filename, bool cache = false);

  /*! Loads image from JPEG file. */
  Ref<Image> loadJPEG(const FileName& fileName);

  /*! Loads image using ImageMagick. */
  Ref<Image> loadMagick(const FileName& fileName);
  
  /*! Loads image from PFM file. */
  Ref<Image> loadPFM(const FileName& fileName);

  /*! Loads image from PNG file. */
//Ref<Image> loadPNG(const FileName& fileName);

  /*! Loads image from PPM file. */
  Ref<Image> loadPPM(const FileName& fileName);

  /*! Loads image from TIFF file. */
//Ref<Image> loadTIFF(const FileName& fileName);
  
  /*! Store image to EXR file. */
  void storeExr(const Ref<Image>& img, const FileName& fileName);
  
  /*! Store image to file. Format is auto detected. */
  void storeImage(const Ref<Image>& image, const FileName& filename);

  /*! Store image to JPEG file. */
  void storeJPEG(const Ref<Image>& img, const FileName& fileName);

  /*! Store image to file using ImageMagick. */
  void storeMagick(const Ref<Image>& img, const FileName& fileName);
  
  /*! Store image to PFM file. */
  void storePFM(const Ref<Image>& img, const FileName& fileName);
  
  /*! Store image to PNG file. */
//void storePNG(const Ref<Image>& img, const FileName& fileName);

  /*! Store image to PPM file. */
  void storePPM(const Ref<Image>& img, const FileName& fileName);
  
  /*! Store image to TGA file. */
  void storeTga(const Ref<Image>& img, const FileName& fileName);

  /*! Store image to TIFF file. */
//void storeTIFF(const Ref<Image>& img, const FileName& fileName);

}
