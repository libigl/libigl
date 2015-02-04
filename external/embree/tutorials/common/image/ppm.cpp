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

#include "image/image.h"

#include <iostream>
#include <cstring>
#include <cstdio>

namespace embree
{
  /*! read a single comment line starting with #, or read a space */
  static bool readCommentLine(FILE* file)
  {
    int c = fgetc(file);
    if (isspace(c)) return true;
    if (c != '#') {
      ungetc(c, file);
      return false;
    }
    char line[1024];
    if(fgets(line, sizeof(line), file) == NULL) 
      THROW_RUNTIME_ERROR("Error reading PPM file!");
    return true;
  }

  /*! read PPM file from disk */
  Ref<Image> loadPPM(const FileName& fileName)
  {
    /* open PPM file */
    FILE* file = fopen(fileName.c_str(), "rb");
    if (!file) THROW_RUNTIME_ERROR("cannot open " + fileName.str());

    /* read file type */
    char type[8];
    if (fscanf(file, "%7s", type) != 1)
      THROW_RUNTIME_ERROR("Error reading " + fileName.str());
    
    /* skip comment lines */
    while (readCommentLine(file)) {};

    /* read width, height, and maximal color value */
    int width, height, maxColor;
    if (fscanf(file, "%i %i %i", &width, &height, &maxColor) != 3)
      THROW_RUNTIME_ERROR("Error reading " + fileName.str());
    float rcpMaxColor = 1.0f/float(maxColor);

    /* get return or space */
    fgetc(file);

    /* create image and fill with data */
    Ref<Image> img = new Image4c(width,height,fileName);

    /* image in text format */
    if (!strcmp(type, "P3"))
    {
      int r, g, b;
      for (ssize_t y=0; y<height; y++) {
        for (ssize_t x=0; x<width; x++) {
          if (fscanf(file, "%i %i %i", &r, &g, &b) != 3)
            THROW_RUNTIME_ERROR("Error reading " + fileName.str());
          img->set(x,y,Color4(float(r)*rcpMaxColor,float(g)*rcpMaxColor,float(b)*rcpMaxColor,1.0f));
        }
      }
    }

    /* image in binary format 8 bit */
    else if (!strcmp(type, "P6") && maxColor <= 255)
    {
      unsigned char rgb[3];
      for (ssize_t y=0; y<height; y++) {
        for (ssize_t x=0; x<width; x++) {
          if (fread(rgb,sizeof(rgb),1,file) != 1)
            THROW_RUNTIME_ERROR("Error reading " + fileName.str());
          img->set(x,y,Color4(float(rgb[0])*rcpMaxColor,float(rgb[1])*rcpMaxColor,float(rgb[2])*rcpMaxColor,1.0f));
        }
      }
    }

    /* image in binary format 16 bit */
    else if (!strcmp(type, "P6"))
    {
      unsigned short rgb[3];
      for (ssize_t y=0; y<height; y++) {
        for (ssize_t x=0; x<width; x++) {
          if (fread(rgb,sizeof(rgb),1,file) != 1)
            THROW_RUNTIME_ERROR("Error reading " + fileName.str());
          img->set(x,y,Color4(float(rgb[0])*rcpMaxColor,float(rgb[1])*rcpMaxColor,float(rgb[2])*rcpMaxColor,1.0f));
        }
      }
    }

    /* invalid magic value */
    else {
      fclose(file);
      THROW_RUNTIME_ERROR("Invalid magic value in PPM file");
    }

    fclose(file);
    return img;
  }

  /*! store PPM file to disk */
  void storePPM(const Ref<Image>& img, const FileName& fileName)
  {
    FILE* file = fopen(fileName.c_str(), "wb");
    if (!file) THROW_RUNTIME_ERROR("cannot open file " + fileName.str());
    fprintf(file,"P6\n%i %i\n255\n", int(img->width), int(img->height));

    for (size_t y=0; y<img->height; y++) {
      for (size_t x=0; x<img->width; x++) {
        const Color4 c = img->get(x,y);
        fputc((unsigned char)(clamp(c.r)*255.0f), file);
        fputc((unsigned char)(clamp(c.g)*255.0f), file);
        fputc((unsigned char)(clamp(c.b)*255.0f), file);
      }
    }
    fclose(file);
  }
}
