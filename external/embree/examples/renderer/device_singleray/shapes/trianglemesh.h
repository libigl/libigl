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

#ifndef __EMBREE_TRIANGLE_MESH_H__
#define __EMBREE_TRIANGLE_MESH_H__

#include "trianglemesh_normals.h"
#include "trianglemesh_full.h"

namespace embree
{
  class TriangleMesh
  {
  public:

    static Ref<Shape> create (const Parms& parms) 
    {
      bool hasPositions = parms.getData("positions");
      bool hasMotions   = parms.getData("motions");
      bool hasNormals   = parms.getData("normals");
      bool hasTangents  = parms.getData("tangent_x") | parms.getData("tangent_y");
      bool hasTexCoords = parms.getData("texcoords") | parms.getData("texcoords0");

      if (hasPositions && !hasMotions && hasNormals && !hasTangents && !hasTexCoords)
        return new TriangleMeshWithNormals(parms);
      else
        return new TriangleMeshFull(parms);
    }
  };
}

#endif
