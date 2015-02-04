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

#include <cstdio>

namespace embree
{
  inline void fwrite_uchar (unsigned char  v, FILE* file) { fwrite(&v, sizeof(v), 1, file); }
  inline void fwrite_ushort(unsigned short v, FILE* file) { fwrite(&v, sizeof(v), 1, file); }

  void storeTga(const Ref<Image>& img, const FileName& fileName)
  {
    FILE* file = fopen(fileName.c_str(), "wb");
    if (!file) THROW_RUNTIME_ERROR("error opening file " + fileName.str());

    fwrite_uchar(0x00, file);
    fwrite_uchar(0x00, file);
    fwrite_uchar(0x02, file);
    fwrite_ushort(0x0000, file);
    fwrite_ushort(0x0000, file);
    fwrite_ushort(0x0000, file);
    fwrite_ushort(0x0000, file);
    fwrite_uchar(0x00, file);
    fwrite_ushort((unsigned short)img->width , file);
    fwrite_ushort((unsigned short)img->height, file);
    fwrite_uchar(0x18, file);
    fwrite_uchar(0x20, file);

    for (size_t y=0; y<img->height; y++) {
      for (size_t x=0; x<img->width; x++) {
        Color c = img->get(x,y);
        fwrite_uchar((unsigned char)(clamp(c.b)*255.0f), file);
        fwrite_uchar((unsigned char)(clamp(c.g)*255.0f), file);
        fwrite_uchar((unsigned char)(clamp(c.r)*255.0f), file);
      }
    }
    fclose(file);
  }
}
