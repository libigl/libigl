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

#include "triangle4.h"
#if defined(__TARGET_AVX__)
#include "triangle8.h"
#endif
#include "common/scene.h"

namespace embree
{
  Triangle4Type Triangle4Type::type;
  TriangleMeshTriangle4 TriangleMeshTriangle4::type;

  Triangle4Type::Triangle4Type () 
  : PrimitiveType("triangle4",sizeof(Triangle4),4,false,1) {} 

#if defined(__TARGET_AVX__)
  Triangle8Type Triangle8Type::type;
  TriangleMeshTriangle8 TriangleMeshTriangle8::type;

  Triangle8Type::Triangle8Type () 
    : PrimitiveType("triangle8",2*sizeof(Triangle4),8,false,1) {}
#endif
  
  size_t Triangle4Type::blocks(size_t x) const {
    return (x+3)/4;
  }
  
  size_t Triangle4Type::size(const char* This) const {
    return ((Triangle4*)This)->size();
  }
  
  size_t Triangle4Type::hash(const char* This, size_t num) const 
  {
    size_t hash = 0;
    for (size_t i=0; i<num; i++)
      hash += (i+1)*((Triangle4*)This)[i].hash();
    return hash;
  }
  
  BBox3fa TriangleMeshTriangle4::update(char* prim_i, size_t num, void* geom) const 
  {
    BBox3fa bounds = empty;
    TriangleMesh* mesh = (TriangleMesh*) geom;
    Triangle4* prim = (Triangle4*)prim_i;

    if (num == -1)
    {
      while (true)
      {
	ssei vgeomID = -1, vprimID = -1, vmask = -1;
	sse3f v0 = zero, v1 = zero, v2 = zero;
	
	for (size_t i=0; i<4; i++)
	{
	  if (prim->geomID<1>(i) == -1) break;
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
	new (prim) Triangle4(v0,v1,v2,vgeomID,vprimID,vmask,last);
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
	  if (prim->geomID<0>(i) == -1) break;
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
	new (prim) Triangle4(v0,v1,v2,vgeomID,vprimID,vmask,false);
      }
    }
    return bounds; 
  }
}
