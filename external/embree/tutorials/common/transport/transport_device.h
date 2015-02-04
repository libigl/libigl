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

#include "math/vec3.h"

namespace embree
{
  struct Scene
  {
    void* materials;  //!< material list
    void* positions;    //!< vertex position array
    void* normals;       //!< vertex normal array
    void* texcoords;     //!< vertex texcoord array
    void* triangles;  //!< list of triangles
    int numMaterials;
    int numVertices;
    int numTriangles;
    embree::Vec3f pointLightPosition;
    embree::Vec3f pointLightIntensity;
    embree::Vec3f ambientLightIntensity;
  };
  
  extern "C" void device_init(const char* cfg);
  extern "C" void device_key_pressed (int32 key);

  extern "C" void device_set_scene(Scene* scene);
  extern "C" void device_resize(int32 width, int32 height);
  extern "C" bool device_pick(const float x, const float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p, Vec3fa& hitPos);

  extern "C" void device_render(int* pixels, const int width, const int height,
                                const float time, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p);
  extern "C" void device_cleanup();
}
