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

#include "triangle4i.h"
#include "common/scene.h"

namespace embree
{
  Triangle4iType Triangle4iType::type;
  TriangleMeshTriangle4i TriangleMeshTriangle4i::type;

  Triangle4iType::Triangle4iType () 
    : PrimitiveType("triangle4i",sizeof(Triangle4i),4,true,1) {} 
  
  size_t Triangle4iType::blocks(size_t x) const {
    return (x+3)/4;
  }
  
  size_t Triangle4iType::size(const char* This) const {
    return ((Triangle4i*)This)->size();
  }

  BBox3fa TriangleMeshTriangle4i::update(char* prim_i, size_t num, void* geom) const 
  {
    BBox3fa bounds = empty;
    TriangleMesh* mesh = (TriangleMesh*) geom;
    Triangle4i* prim = (Triangle4i*) prim_i;

    if (num == -1)
    {
      while (true)
      {
	ssei vgeomID = -1, vprimID = -1, vmask = -1;
	sse3f v0 = zero, v1 = zero, v2 = zero;
	
	for (size_t i=0; i<4; i++)
	{
	  if (prim->primID<1>(i) == -1) break;
	  const unsigned geomID = prim->geomID<1>(i);
	  const unsigned primID = prim->primID<1>(i);
	  const TriangleMesh::Triangle& tri = mesh->triangle(primID);
	  const Vec3fa p0 = mesh->vertex(tri.v[0]);
	  const Vec3fa p1 = mesh->vertex(tri.v[1]);
	  const Vec3fa p2 = mesh->vertex(tri.v[2]);
	  bounds.extend(merge(BBox3fa(p0),BBox3fa(p1),BBox3fa(p2)));
	}
	bool last = prim->last();
	if (last) break;
	prim++;
      }
    }
    else
    {
      for (size_t j=0; j<num; j++, prim++)
      {
	ssei vgeomID = -1, vprimID = -1, vmask = -1;
	sse3f v0 = zero, v1 = zero, v2 = zero;
	
	for (size_t i=0; i<4; i++)
	{
	  if (prim->primID<0>(i) == -1) break;
	  const unsigned geomID = prim->geomID<0>(i);
	  const unsigned primID = prim->primID<0>(i);
	  const TriangleMesh::Triangle& tri = mesh->triangle(primID);
	  const Vec3fa p0 = mesh->vertex(tri.v[0]);
	  const Vec3fa p1 = mesh->vertex(tri.v[1]);
	  const Vec3fa p2 = mesh->vertex(tri.v[2]);
	  bounds.extend(merge(BBox3fa(p0),BBox3fa(p1),BBox3fa(p2)));
	}
      }
    }
    return bounds; 
  }
}
