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

#ifndef __EMBREE_ISPC_TRIANGLE_H__
#define __EMBREE_ISPC_TRIANGLE_H__

#include "api/parms.h"
#include "trianglemesh_ispc.h"

namespace embree
{
  class Triangle
  {
  public:

    /*! Construction from parameter container. */
    static void* create (const Parms& parms) 
    {
      Vector3f v0 = parms.getVector3f("v0");
      Vector3f v1 = parms.getVector3f("v1");
      Vector3f v2 = parms.getVector3f("v2");
      Vector3f Ng = normalize(cross(v2-v0,v1-v0));

      int numVertices = 3;
      Vec3fa* position = new Vec3fa[3];
      position[0] = v0;
      position[1] = v1;
      position[2] = v2;
      Vec3fa* motion = NULL;
      Vec3fa* normal = new Vec3fa[3];
      normal[0] = Ng;
      normal[1] = Ng;
      normal[2] = Ng;
      Vec2f* texcoord = new Vec2f[3];
      texcoord[0] = Vec2f(0.0f,0.0f);
      texcoord[1] = Vec2f(0.0f,1.0f);
      texcoord[2] = Vec2f(1.0f,0.0f);
      int numTriangles = 1;
      Vec4i* triangles = new Vec4i[1];
      triangles[0] = Vec4i(0,1,2,0);
    
      return  ispc::TriangleMesh__new((ispc::vec3fa*)position,
                                      (ispc::vec3fa*)motion,
                                      (ispc::vec3fa*)normal,
                                      (ispc::vec2f*)texcoord,
                                      numVertices,
                                      (ispc::vec4i*)triangles,
                                      numTriangles);
    }
  };
}

#endif
