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

#ifndef __EMBREE_DISK_H__
#define __EMBREE_DISK_H__

#include "../shapes/trianglemesh_full.h"

namespace embree
{
  /*! Implements a disk shape. */
  class Disk : public TriangleMeshFull
  {
  public:

    /*! Construction from acceleration structure type. */
    Disk (const AccelType& ty)
      : TriangleMeshFull(ty) {}

    /*! Construction from parameter container. */
    Disk (const Parms& parms) 
      : TriangleMeshFull(AccelType(parms))
    {
      P    = parms.getVector3f("P");
      h    = parms.getFloat("h");
      r    = parms.getFloat("r");
      numTriangles = parms.getInt("numTriangles");
      triangulate();
    }

  private:

    void triangulate()
    {
      float rcpNumTriangles = rcp(float(numTriangles));
      for (size_t phi=0; phi<numTriangles; phi++)
      {
        Vector3f d(sinf(phi*2.0f*float(pi)*rcpNumTriangles),cosf(phi*2.0f*float(pi)*rcpNumTriangles),0.0f);
        position.push_back(P+r*d);
        normal.push_back(Vector3f(0.0f,0.0f,1.0f));
        texcoord.push_back(Vec2f(0.0f,0.0f));
      }
      position.push_back(P+Vector3f(0,0,h));
      for (size_t phi=0; phi<numTriangles; phi++) {
        size_t p0 = numTriangles;
        size_t p1 = (phi+0)%numTriangles;
        size_t p2 = (phi+1)%numTriangles;
        switch (phi%3) {
        case 0: triangles.push_back(Vector3i((int)p0,(int)p2,(int)p1)); break;
        case 1: triangles.push_back(Vector3i((int)p1,(int)p0,(int)p2)); break;
        case 2: triangles.push_back(Vector3i((int)p2,(int)p1,(int)p0)); break;
        }
      }
    }

  public:
    Vector3f P;           //!< center of disk
    float h;           //!< height of the cone
    float r;           //!< radius of disk
    size_t numTriangles;   //!< triangulation amount
  };
}

#endif
