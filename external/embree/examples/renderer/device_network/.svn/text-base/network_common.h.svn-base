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

#ifndef __EMBREE_NETWORK_COMMON_H__
#define __EMBREE_NETWORK_COMMON_H__

#include "default.h"
#include "sys/network.h"

namespace embree
{
  /*! magick code of network packets */
  static const int magick = 0x32657845;

  /*! enum type containing all network commands */
  enum {
    EMBREE_NEW_CAMERA = 0,
    EMBREE_NEW_DATA = 1,
    EMBREE_NEW_DATA_FROM_FILE = 2,
    EMBREE_NEW_IMAGE = 3,
    EMBREE_NEW_IMAGE_FROM_FILE = 4,
    EMBREE_NEW_TEXTURE = 5,
    EMBREE_NEW_MATERIAL = 6,
    EMBREE_NEW_SHAPE = 7,
    EMBREE_NEW_LIGHT = 8,
    EMBREE_NEW_SHAPE_PRIMITIVE = 9,
    EMBREE_NEW_LIGHT_PRIMITIVE = 10,
    EMBREE_TRANSFORM_PRIMITIVE = 44,
    EMBREE_NEW_SCENE = 11,
    EMBREE_NEW_TONEMAPPER = 12,
    EMBREE_NEW_RENDERER = 13,
    EMBREE_NEW_FRAMEBUFFER = 14,
    EMBREE_DECREF = 15,
    EMBREE_INCREF = 16,
    EMBREE_SET_BOOL1 = 17,
    EMBREE_SET_BOOL2 = 18,
    EMBREE_SET_BOOL3 = 19,
    EMBREE_SET_BOOL4 = 20,
    EMBREE_SET_INT1 = 21,
    EMBREE_SET_INT2 = 22,
    EMBREE_SET_INT3 = 23,
    EMBREE_SET_INT4 = 24,
    EMBREE_SET_FLOAT1 = 25,
    EMBREE_SET_FLOAT2 = 26,
    EMBREE_SET_FLOAT3 = 27,
    EMBREE_SET_FLOAT4 = 28,
    EMBREE_SET_ARRAY = 29,
    EMBREE_SET_STRING = 30,
    EMBREE_SET_IMAGE = 31,
    EMBREE_SET_TEXTURE = 32,
    EMBREE_SET_TRANSFORM = 33,
    EMBREE_SET_SCENE_PRIMITIVE = 50,
    EMBREE_CLEAR = 34,
    EMBREE_COMMIT = 35,
    EMBREE_RENDER_FRAME = 36,
    EMBREE_SWAP_BUFFERS = 37,
    EMBREE_PICK = 38,
    EMBREE_PICK_RESULT = 39,
    EMBREE_FRAME_DATA_RGB_FLOAT = 40,
    EMBREE_FRAME_DATA_RGB8 = 41,
    EMBREE_FRAME_DATA_RGBE8 = 42,
    EMBREE_FRAME_DATA_NATIVE = 49,
    EMBREE_FRAME_DATA_DXT1 = 43,
    EMBREE_FRAME_DATA_JPEG = 45,
    EMBREE_RENDER_TIME = 46
  };

  /*! RGBE encoding */
  __forceinline int encodeRGBE8(const Vector3f& v) {
    float largest = max(v.x,v.y,v.z);
    int bound = cast_f2i(largest) & 0x7F800000;
    int exp = (bound >> 23) - 7;
    int nexp = (-exp+254) << 23;
    float scale = cast_i2f(nexp);
    Vector3f vs = scale*v;
    return (int(vs.x) << 0) | (int(vs.y) << 8) | (int(vs.z) << 16) | (int(exp) << 24); 
  }

  /*! RGBE decoding */
  __forceinline Vector3f decodeRGBE8(int rgbe) {
    unsigned char r = (unsigned char)(rgbe >>  0);
    unsigned char g = (unsigned char)(rgbe >>  8);
    unsigned char b = (unsigned char)(rgbe >> 16);
    unsigned char e = (unsigned char)(rgbe >> 24);
    Vector3f rgb = Vector3f(float(r),float(g),float(b));
    int exp = int(e) << 23;
    float mscale = cast_i2f(exp);
    return rgb * mscale;
  }
}

#endif
