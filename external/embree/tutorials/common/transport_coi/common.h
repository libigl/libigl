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

namespace embree
{
  struct InitData 
  {
    char cfg[1024];
  };

  struct KeyPressedData
  {
    int32 key;
  };

  struct CreateSceneData 
  {
    int numMaterials;
    int numMeshes;
    //Vec3f pointLightPosition;
    //Vec3f pointLightIntensity;
    //Vec3f ambientLightIntensity;
  };

  struct CreateMeshData
  {
    int numVertices;
    int numTriangles;
    Vec3f dir;
    Vec3f offset;
  };

  struct ResizeData {
    int32 width, height;
  };

  struct PickDataSend {
    float x,y;
    Vec3f vx,vy,vz,p;
  };

  struct PickDataReceive {
    Vec3f pos;
    bool hit;
  };

  struct RenderData {
    float time;
    Vec3f vx,vy,vz,p;
    int width;
    int height;
  };
}
