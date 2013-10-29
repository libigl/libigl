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

#ifndef __EMBREE_TUTORIAL03_OBJ_LOADER_H__
#define __EMBREE_TUTORIAL03_OBJ_LOADER_H__

#include "sys/platform.h"
#include "sys/filename.h"
#include "math/vec2.h"
#include "math/vec3.h"

#include <vector>

namespace embree
{
  /*! mesh representing the OBJ file */
  struct OBJMesh 
  {
    /*! OBJ Triangle */
    struct Triangle 
    {
    public:
      Triangle (int v0, int v1, int v2, int triangleID, int materialID) 
      : v0(v0), v1(v1), v2(v2), triangleID(triangleID), materialID(materialID) {}

    public:
      int v0, v1, v2, triangleID, materialID;
    };
    
    /*! OBJ material */
    struct Material
    {
    public:
      Material ()
      : illum(0), d(1.f), Ns(1.f), Ni(1.f), Ka(0.f), Kd(1.f), Ks(0.f), Tf(1.f) {};
      
    public:
      int illum;             /*< illumination model */
      
      float d;               /*< dissolve factor, 1=opaque, 0=transparent */
      float Ns;              /*< specular exponent */
      float Ni;              /*< optical density for the surface (index of refraction) */
      
      Vec3f Ka;              /*< ambient reflectivity */
      Vec3f Kd;              /*< diffuse reflectivity */
      Vec3f Ks;              /*< specular reflectivity */
      Vec3f Tf;              /*< transmission filter */
    };

  public:
    std::vector<Material> material;   //!< material list
    std::vector<Vec3fa> v;	      //!< vertex position array
    std::vector<Vec3f> vn;            //!< vertex normal array
    std::vector<Vec2f> vt;            //!< vertex texcoord array
    std::vector<Triangle> triangles;  //!< list of triangles
  };
  
  /*! read from disk */
  void loadOBJ(const FileName& fileName, OBJMesh& mesh);
}

#endif
