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

#ifdef USE_IMAGEMAGICK

#include "image/image.h"
#include <iostream>

/*! include Image Magick headers */
#include <Magick++.h>
using namespace Magick;

namespace embree
{
  Ref<Image> loadMagick(const FileName& fileName)
  {
    Magick::Image image(fileName.c_str());
    Image* out = new Image4c(image.columns(),image.rows(),fileName);
    float rcpMaxRGB = 1.0f/float(MaxRGB);
    Magick::Pixels pixel_cache(image);
    Magick::PixelPacket* pixels = pixel_cache.get(0,0,out->width,out->height);
    
    for (size_t y=0; y<out->height; y++) {
      for (size_t x=0; x<out->width; x++) {
        Color4 c;
        c.r = float(pixels[y*out->width+x].red    )*rcpMaxRGB;
        c.g = float(pixels[y*out->width+x].green  )*rcpMaxRGB;
        c.b = float(pixels[y*out->width+x].blue   )*rcpMaxRGB;
        c.a = float(pixels[y*out->width+x].opacity)*rcpMaxRGB;
        out->set(x,y,c);
      }
    }

    return out;
  }

  void storeMagick(const Ref<Image>& img, const FileName& fileName)
  {
    Magick::Image image(Magick::Geometry(img->width,img->height),
                        Magick::ColorRGB(0,0,0));
    image.modifyImage();

    Magick::Pixels pixel_cache(image);
    Magick::PixelPacket* pixels = pixel_cache.get(0,0,img->width,img->height);
    for (size_t y=0; y<img->height; y++) {
      for (size_t x=0; x<img->width; x++) {
        Color4 c = img->get(x,y);
	pixels[y*img->width+x].red     = Magick::Quantum(clamp(c.r)*MaxRGB);
	pixels[y*img->width+x].green   = Magick::Quantum(clamp(c.g)*MaxRGB);
	pixels[y*img->width+x].blue    = Magick::Quantum(clamp(c.b)*MaxRGB);
	pixels[y*img->width+x].opacity = Magick::Quantum(clamp(c.a)*MaxRGB);
      }
    }
    pixel_cache.sync();
    image.write(fileName.c_str());
  }
}

#endif // USE_IMAGEMAGICK

