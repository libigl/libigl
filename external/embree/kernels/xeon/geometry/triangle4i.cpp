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

#include "triangle4i.h"
#include "common/scene.h"

namespace embree
{
  Triangle4iType Triangle4iType::type;

  Triangle4iType::Triangle4iType () 
  : PrimitiveType("triangle4i",sizeof(Triangle4i),4,true,1) {} 
  
  size_t Triangle4iType::blocks(size_t x) const {
    return (x+3)/4;
  }
  
  size_t Triangle4iType::size(const char* This) const {
    return ((Triangle4i*)This)->size();
  }
  
  void Triangle4iType::pack(char* This, atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, void* geom) const 
  {
    Scene* scene = (Scene*) geom;
    
    ssei geomID = -1, primID = -1;
    Vec3f* v0[4] = { NULL, NULL, NULL, NULL };
    ssei v1 = zero, v2 = zero;
    PrimRef& prim = *prims;
    
    for (size_t i=0; i<4; i++)
    {
      const TriangleMeshScene::TriangleMesh* mesh = scene->getTriangleMesh(prim.geomID());
      const TriangleMeshScene::TriangleMesh::Triangle& tri = mesh->triangle(prim.primID());
      if (prims) {
        geomID[i] = prim.geomID();
        primID[i] = prim.primID();
        v0[i] = (Vec3f*) &mesh->vertex(tri.v[0]); 
        v1[i] = (int*)&mesh->vertex(tri.v[1])-(int*)v0[i]; 
        v2[i] = (int*)&mesh->vertex(tri.v[2])-(int*)v0[i]; 
        prims++;
      } else {
        assert(i);
        geomID[i] = -1;
        primID[i] = -1;
        v0[i] = v0[i-1];
        v1[i] = 0; 
        v2[i] = 0;
      }
      if (prims) prim = *prims;
    }
    
    new (This) Triangle4i(v0,v1,v2,geomID,primID);
  }
}
