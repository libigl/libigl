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

#ifndef __EMBREE_COI_COMMON_H__
#define __EMBREE_COI_COMMON_H__

#include "sys/platform.h"

namespace embree
{
  struct parmsNewCamera {
    int id;
    char type[64];
  };

  struct parmsNewData {
    int id;
    size_t bytes;
  };

  struct parmsNewDataStart {
    int id;
    size_t bytes;
  };

  struct parmsNewDataSet {
    int id;
    size_t offset;
    size_t bytes;
  };

  struct parmsNewImage {
    int id;
    size_t width;
    size_t height;
    char type[64];
  };

  struct parmsNewTexture {
    int id;
    char type[64];
  };

  struct parmsNewMaterial {
    int id;
    char type[64];
  };

  struct parmsNewShape {
    int id;
    char type[64];
  };

  struct parmsNewLight {
    int id;
    char type[64];
  };

  struct parmsNewShapePrimitive {
    int id;
    int shape;
    int material;
    float transform[12];
  };

  struct parmsNewLightPrimitive {
    int id;
    int light;
    int material;
    float transform[12];
  };

  struct parmsTransformPrimitive {
    int id;
    int primitive;
    float transform[12];
  };

  struct parmsNewScene {
    int id;
    char type[64];
  };

  struct parmsSetPrimitive {
    int scene;
    int slot;
    int prim;
  };

  struct parmsNewToneMapper {
    int id;
    char type[64];
  };

  struct parmsNewRenderer {
    int id;
    char type[64];
  };

  struct parmsNewFrameBuffer {
    int id;
    size_t width;
    size_t height;
    size_t depth;
    char type[64];
  };

  struct parmsSwapBuffers {
    int framebuffer;
  };

  struct parmsIncRef {
    int handle;
  };

  struct parmsDecRef {
    int handle;
  };

  struct parmsSetBoolN {
    int handle;
    bool x;
    bool y;
    bool z;
    bool w;
    char property[64];
  };

  struct parmsSetIntN {
    int handle;
    int x;
    int y;
    int z;
    int w;
    char property[64];
  };

  struct parmsSetFloatN {
    int handle;
    float x;
    float y;
    float z;
    float w;
    char property[64];
  };

  struct parmsSetArray {
    int handle;
    char property[64];
    char type[64];
    int data;
    size_t size;
    size_t stride;
    size_t ofs;
  };

  struct parmsSetString {
    int handle;
    char property[64];
    char str[64];
  };

  struct parmsSetImage {
    int handle;
    int image;
    char property[64];
  };

  struct parmsSetTexture {
    int handle;
    int texture;
    char property[64];
  };

  struct parmsSetTransform {
    int handle;
    float transform[12];
    char property[64];
  };

  struct parmsClear {
    int handle;
  };

  struct parmsCommit {
    int handle;
  };

  struct parmsRenderFrame {
    int renderer;
    int camera;
    int scene;
    int toneMapper;
    int frameBuffer;
    int accumulate;
  };

  struct parmsPick {
    int camera;
    float x,y;
    int scene;
  };

  struct returnPick {
    bool hit;
    float x,y,z;
  };
}

#endif
