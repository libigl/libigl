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

#include "triangle8.h"
#include "common/scene.h"

namespace embree
{
  size_t Triangle8Type::blocks(size_t x) const {
    return (x+7)/8;
  }
  
  size_t Triangle8Type::size(const char* This) const {
    return ((Triangle8*)This)->size();
  }
  
  BBox3fa TriangleMeshTriangle8::update(char* prim_i, size_t num, void* geom) const 
  {
    BBox3fa bounds = empty;
    TriangleMesh* mesh = (TriangleMesh*) geom;
    Triangle8* prim = (Triangle8*) prim_i;

    if (num == -1)
    {
      while (true)
      {
	avxi vgeomID = -1, vprimID = -1, vmask = -1;
	avx3f v0 = zero, v1 = zero, v2 = zero;
	
	for (size_t i=0; i<8; i++)
	{
	  if (prim->primID<1>(i) == -1) break;
	  const unsigned geomID = prim->geomID<1>(i);
	  const unsigned primID = prim->primID<1>(i);
	  const TriangleMesh::Triangle& tri = mesh->triangle(primID);
	  const Vec3fa p0 = mesh->vertex(tri.v[0]);
	  const Vec3fa p1 = mesh->vertex(tri.v[1]);
	  const Vec3fa p2 = mesh->vertex(tri.v[2]);
	  bounds.extend(merge(BBox3fa(p0),BBox3fa(p1),BBox3fa(p2)));
	  vgeomID [i] = geomID;
	  vprimID [i] = primID;
	  vmask   [i] = mesh->mask;
	  v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
	  v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
	  v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
	}
	bool last = prim->last();
	new (prim) Triangle8(v0,v1,v2,vgeomID,vprimID,vmask,last);
	if (last) break;
	prim++;
      }
    }
    else
    {
      for (size_t j=0; j<num; j++, prim++)
      {
	avxi vgeomID = -1, vprimID = -1, vmask = -1;
	avx3f v0 = zero, v1 = zero, v2 = zero;
	
	for (size_t i=0; i<8; i++)
	{
	  if (prim->primID<0>(i) == -1) break;
	  const unsigned geomID = prim->geomID<0>(i);
	  const unsigned primID = prim->primID<0>(i);
	  const TriangleMesh::Triangle& tri = mesh->triangle(primID);
	  const Vec3fa p0 = mesh->vertex(tri.v[0]);
	  const Vec3fa p1 = mesh->vertex(tri.v[1]);
	  const Vec3fa p2 = mesh->vertex(tri.v[2]);
	  bounds.extend(merge(BBox3fa(p0),BBox3fa(p1),BBox3fa(p2)));
	  vgeomID [i] = geomID;
	  vprimID [i] = primID;
	  vmask   [i] = mesh->mask;
	  v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
	  v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
	  v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
	}
	new (prim) Triangle8(v0,v1,v2,vgeomID,vprimID,vmask,false);
      }
    }
    return bounds; 
  }
}
