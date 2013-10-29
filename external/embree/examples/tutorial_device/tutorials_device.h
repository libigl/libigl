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

#include "math/vec3.h"

namespace ispc 
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

  extern "C" void init(int verbose);
  extern "C" void set_scene(Scene* scene);
  extern "C" void resize(int32 width, int32 height);
  extern "C" void render(void* ptr, const int width, const int height,
                         const float time, const struct vec3f& vx, const struct vec3f& vy, const struct vec3f& vz, const struct vec3f& p);
  extern "C" void cleanup();
}
