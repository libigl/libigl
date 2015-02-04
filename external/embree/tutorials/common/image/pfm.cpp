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
    if (fgets(line, sizeof(line), file) == NULL)
      THROW_RUNTIME_ERROR("Error reading PPM file!");
    return true;
  }

  /*! read PFM file from disk */
  Ref<Image> loadPFM(const FileName& fileName)
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
    int width, height;
    float maxColor;
    if (fscanf(file, "%i %i %f", &width, &height, &maxColor) != 3)
      THROW_RUNTIME_ERROR("Error reading " + fileName.str());

    /* Check for big endian PFM file */
    if (maxColor > 0.0f) {
      fclose(file);
      THROW_RUNTIME_ERROR("Big endian PFM files not supported");
    }
    float rcpMaxColor = -1.0f/float(maxColor);

    /* get return or space */
    fgetc(file);

    /* create image and fill with data */
    Ref<Image> img = new Image3f(width,height,fileName);

    /* image in binary format 32 bit */
    if (!strcmp(type, "PF"))
    {
      float rgb[3];
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
      THROW_RUNTIME_ERROR("Invalid magic value in PFM file");
    }

    fclose(file);
    return img;
  }

  /*! store PFM file to disk */
  void storePFM(const Ref<Image>& img, const FileName& fileName)
  {
    FILE* file = fopen(fileName.c_str(), "wb");
    if (!file) THROW_RUNTIME_ERROR("cannot open file " + fileName.str());
    fprintf(file,"PF\n%i %i\n%f\n", int(img->width), int(img->height), -1.0f);

    float rgb[3];
    for (size_t y=0; y<img->height; y++) {
      for (size_t x=0; x<img->width; x++) {
        Color4 c = img->get(x,y);
        rgb[0] = c.r; rgb[1] = c.g; rgb[2] = c.b;
        fwrite(rgb,sizeof(rgb),1,file);
      }
    }
    fclose(file);
  }
}
